function [cfg]=bg_SWM(cfg, dat)
% [cfg]=bg_SWM(cfg, dat)
% cfg:  a structure that contains all parameters
% dat:  optional the data on which you want to apply the algorithm
%       (not neccesary if cfg.fname and cfg.varname are present)
%       Dimensions: M x N
%       M: trials
%       N: time points
% 
% Fields in cfg that are required:
% .fitlen: the length of the sliding windows in timestamp units
% 
% Optional fields (they all have a default value):
% CORE PARAMETERS
% .fname:     path and filename of the .mat-file that contains the data on
%             which the algorithms must be applied.
% .varname:   the name of the variable within .fname tha actually contains
%             the data.
% .numIt:     the number of iterations the algorithm will run through.
%             (default = 1e4)
% .guard:     the minimal space between the starting positions of two 
%             succesive windows. Also determines the number of windows.
%             (default = .fitlen/1.5)
% .numWin:    number of windows to fit per trial. Limited by length of
%             data, .fitlen and .guard;
%             (default = max)
% .nPT:       number of parallel temperatures 
%             (default = 1)
% .Tfac:      Temperatures used in the PT sampling algorithm. This overrides
%             .nPT to numel(.Tfac). 
%             (default = logspace(-3,1,.nPT))
% .konstant:  Bolzmann constant used in determining whether two
%             temperatures should be exchanged.
%             (default = 1e3)
% 
% Note: good values for .Tfac and .konstant depend on the scaling of your
% data and the signal to noise level
% 
% OPTIONAL PARAMETERS
% .numclust:  number of shapes (=clusters) to find
%             (default = 1)
% .fs:        Sampling rate of dataset. Needed when using .Fbp or .Fbs
%             below
% .Fbp:       Bandpass frequency range. Usage: [hp lp]
%             When .Fbp is supplied, the function first applies a bandpass
%             filter to the data before performing the SWM
% .Fbs:       Bandstop frequency range. Usage: [lp hp]
%             When .Fbs is supplied, the function first applies a bandstop
%             filter to the data before performing the SWM
%
% FLAGS
% .fullOutput:Flag to determine wheter the function should output the full
%             trajectory of the cost function, or just a undersampled 
%             version of it. Useful for debugging or checking the effect of
%             chaning the temperature parameters
%             (default = 0)
% .debug:     Flag to determine wheter the algorithm should run in debug
%             mode. In debug mode, algorithm will only initilize one
%             temperature and then overwrite all the others. This speeds up
%             initilization of the algorithm and is therefore useful during
%             debugging.


% load data from file instead of from function input
if nargin<2
  dum=load(cfg.fname, cfg.varname);
  eval(['dat=dum.' cfg.varname ';'])
else
  cfg.fname='function input';
  cfg.varname='N/A';
end


if iscolumn(dat);
  dat=dat.';
end

if isfield(cfg,'Fbp')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. Bandpass filter is not possible without cfg.fs')
  end
  Fbp=cfg.Fbp;
  dat=ft_preproc_bandpassfilter(dat, fs, Fbp);
end

if isfield(cfg,'Fbs')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. Bandstop filter is not possible without cfg.fs')
  end
  Fbs=cfg.Fbs;
  dat=ft_preproc_bandstopfilter(dat, fs, Fbs);
end

% dat : mxn : m trials, n time points
% 2*guard should be bigger than fitlen
if isfield(cfg,'fitlen')
  fitlen=cfg.fitlen;
else
  error('fitlen should be defined')
end

if isfield(cfg,'guard')
  guard=cfg.guard;
else
  guard=round(fitlen/1.5);
  cfg.guard=guard;
end


if isfield(cfg,'numIt')
  numIt=cfg.numIt;
else
  numIt=1e4;
  cfg.numIt=numIt;
end

% add CoM weighting (not recommended; biases towards symmetry)
if isfield(cfg,'ratio')
  ratio=cfg.ratio;
else
  ratio=0;
  cfg.ratio=ratio;
end

if isfield(cfg,'konstant')
  konstant=cfg.konstant;
else
  konstant=1e3;
  cfg.konstant=konstant;
end

if isfield(cfg,'debug')
  debug=cfg.debug;
else
  debug=false;
end
  

if isfield(cfg,'nPT')
  nPT=cfg.nPT;
else
  nPT=1;
end

if isfield(cfg,'Tfac')
  Tfac=cfg.Tfac;
  nPT=numel(Tfac);
  cfg.nPT=nPT;
else
  Tfac=logspace(-3,1,nPT);
end

if isfield(cfg,'numWin')
  numWin=cfg.numWin;
else
  numWin=0;
end

if isfield(cfg,'loc')
  loc=cfg.loc;
  if numel(loc)~= nPT || size(loc{1},1)~=size(dat,1)
    error('location is not correct')
  end
else
  loc=cell(1,nPT);
  for n=1:nPT
    if ~debug || n<2
      if numWin
      [loc{n}, numWin]=initloc(guard,fitlen,dat,numWin);
      else
        [loc{n}, numWin]=initloc(guard,fitlen,dat);
      end
    else
      loc{n}=loc{end};
    end
  end
  tloc=loc{end};
  cfg.numWin=numWin;
end

if isfield(cfg,'bestloc')
  loc{end}=cfg.bestloc;
end

if isfield(cfg,'numclust')
  numclust=cfg.numclust;
else
  cfg.numclust=1;
  numclust=1;
end

if isfield(cfg,'clust')
  initclust=0;
  clust=cfg.clust;
else
  initclust=1;
end

if isfield(cfg,'verbose')
  verbose=cfg.verbose;
else
  verbose=true;
end

if isfield(cfg,'fullOutput')
  fullOutput=cfg.fullOutput;
else
  fullOutput=false;
  cfg.fullOutput=fullOutput;
end




% find the number of sample vectors/templates
numtrl=size(dat,1);
numtemplates=numWin*numtrl;
cfg.numtemplates=numtemplates;

% construct sample vectors
z_i=cell(1,nPT);
z_i(:)={nan(numtemplates,fitlen)};
icomcost=nan(numtemplates,nPT);

cm=icomcost;
comcost=nan(1,nPT);
if verbose
  fprintf('\nInitializing sample vectors...')
end
reverseStr='';

for T=1:nPT
  if T<2 || ~debug
  for n=1:numtemplates
    [trldum idxdum]=ind2sub([numtrl,numWin],n);
    s_i=dat(trldum,loc{T}(trldum, idxdum):loc{T}(trldum, idxdum)+fitlen-1);
    z_i{T}(n,:)=(s_i(:)-nanmean(s_i(:)))/nanstd(s_i(:),1);
    
    msg=sprintf('\n Temperature %d \n Template %d/%d', [T n numtemplates]);
    if verbose
      fprintf([reverseStr, msg]);
    end
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
     if ratio>0
      [icomcost(n,T) cm(n,T)]=com(z_i{T}(n,:),ratio);
    else
      icomcost(n,T)=0;
      cm(n,T)=nan;
    end
    
  end
  else
    z_i{T}=z_i{1};
    icomcost(:,T)=icomcost(:,1);
    cm(:,T)=cm(n,1);
  end
  comcost(T)=sum(icomcost(:,T));
end



if initclust
  clust=cell(numclust,nPT);
  clustID=nan(numtemplates,T);
  tempperclust=round(numtemplates/numclust);
  for T=1:nPT
    % initialize clusters
    locidx=randperm(numtemplates);
    
    for clustidx=1:numclust
      
      if clustidx<numclust
        locidxdum=locidx(1:tempperclust);
        locidx=locidx(tempperclust+1:end);
      else
        locidxdum=locidx;
      end
      clustID(locidxdum,T)=clustidx;
      [trldum tidxdum]=ind2sub([numtrl, numWin],locidxdum);
      clust{clustidx,T}.linidx=locidxdum;
      clust{clustidx,T}.trl=trldum;
      clust{clustidx,T}.tidx=tidxdum;
      clust{clustidx,T}.numtemplates=numel(locidxdum);
      
      
    end
  end
  
else %only reconstruct clust ID and initial cost
  clustID=nan(numtemplates,nPT);
  for T=1:nPT
    for n=1:numclust;
      
      clustID(clust{n,T}.linidx,T)=n;
      %     clust{n}.z_isum=sum(z_i(clust{n}.linidx,:),1);
    end
  end
end

if isfield(cfg,'bestclust')
  clust(:,nPT)=cfg.bestclust;
  for n=1:numclust
    clustID(clust{n,nPT}.linidx,nPT)=n;
  end
end



% initinalize the cost function
if verbose
  fprintf('\nInitializing cost matrices...')
end
D=zeros(1,nPT);
N=fitlen;
for T=1:nPT
  D(T)=comcost(T);
  for n=1:numclust
    N_c=clust{n,T}.numtemplates;
    z_isum=sum(z_i{T}(clust{n,T}.linidx,:),1);
    clust{n,T}.z_isum=z_isum;
    clust{n,T}.tempcost=(N_c^2/(N_c-1))-(z_isum*z_isum.')/(N*(N_c-1));
    
    
    D(T)=D(T)+clust{n,T}.tempcost;
  end
end




if verbose
  fprintf('\nDone!\n')
end


swapcount=0;
Tchangecount=0;
iterrecalc=1;
icostmean=nan(nPT,2);
icostsd=icostmean;
Ttoken=zeros(1,nPT);
rejcount=zeros(2,nPT);
clustnumel=nan(1,numclust);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Main optimization loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter=1;
mincost=D;
mincostTot=min(D);
totcost=nan(nPT,numIt);
totcost(:,1)=mincost;

cc=zeros(1,nPT);

reverseStr='';
if verbose
  fprintf('\nFinding templates...\n')
end

while iter<numIt %&&  cc<cclim
  iter=iter+1;
  iterrecalc=iterrecalc+1;
  swapcount=swapcount+1;
  Tchangecount=Tchangecount+1;
  rejcount(2,:)=rejcount(2,:)+1;
  
  %determine to do cluster switch or window shift of sample
  pshift=.5;
  if numclust==1 %only one cluster means: never change clusters
    pshift=2;
  end
  
  for T=1:nPT
    lidx=randval(numtemplates);
    
    clustidx=clustID(lidx,T);
    
    shift=rand<pshift;
    
    if shift %shift a window (no cluster swap)
      [trl tidx]=ind2sub([numtrl,numWin],lidx);
      
      ploc=loc{T}(trl,tidx);
      locchange=0;
      loopcount=0;
      while ~locchange && loopcount<11 %find a possible step
        dir=((rand>0.5)-0.5)*2*randval(floor(guard/2));
        nloc=ploc+dir;
        locchange= nloc>0 && loc{T}(trl,tidx+1)-nloc>=guard;
        
        if tidx>1 && locchange %check guard to the left of the sample as well
          locchange=nloc-loc{T}(trl,tidx-1) >= guard;
        end
        
        loopcount=loopcount+1;
      end
      
      
      
      if locchange
        
        N_c=clust{clustidx,T}.numtemplates;
        %       pcost=D(T);
        pcost=clust{clustidx,T}.tempcost;
        pcomcost=icomcost(lidx,T);
        
        s_dum=dat(trl,nloc:nloc+fitlen-1);
        z_dum=(s_dum-mean(s_dum))/std(s_dum,1);
        
        z_sumdum=clust{clustidx,T}.z_isum-z_i{T}(lidx,:)+z_dum;
        ncost=(N_c^2/(N_c-1))-(z_sumdum*z_sumdum.')/(fitlen*(N_c-1));
        
        if ratio>0
          [ncomcost ncm]=com(z_dum,ratio);
        else
          ncomcost=0;
          ncm=nan;
        end
        
        cVal=(pcost-ncost+pcomcost-ncomcost);
      else %could not find a allowed move -> reject
        cVal=-inf;
      end
      
      %accept/reject window shift
      if cVal>0
        reject=0;
      else
        if rand<exp(1/Tfac(T)*cVal)
          reject=0;
        else
          reject=1;
        end
      end
      
      if ~reject
        %update everything
        loc{T}(trl,tidx)=nloc;
        clust{clustidx,T}.z_isum=z_sumdum;
        clust{clustidx,T}.tempcost=ncost;
        D(T)=D(T)-cVal;
        z_i{T}(lidx,:)=z_dum;
        icomcost(lidx,T)=ncomcost;
        comcost(T)=comcost(T)+ncomcost-pcomcost;
        cm(lidx,T)=ncm;
        
        if mincost(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc{T};
          end
          mincost(T)=D(T);
          cc(T)=0;
        else
          cc(T)=cc(T)+1;
        end
      else
        rejcount(1,T)=rejcount(1,T)+1;
      end
      
      
    else %cluster swap
      relidx=clust{clustidx,T}.linidx==lidx;
%       if sum(relidx)~=1
%         disp('help')
%       end
      nclustidx=randval(numclust);
      while nclustidx == clustidx
        nclustidx=randval(numclust);
      end
      pcost=clust{clustidx,T}.tempcost+clust{nclustidx,T}.tempcost;
      
      %     z_sumdum=[sum(z_i(clust{pclustidx}.linidx(~sel),:),1); sum(z_i([clust{clustidx}.linidx lidx],:),1)];
      z_sumdum=[clust{clustidx,T}.z_isum-z_i{T}(lidx,:); clust{nclustidx,T}.z_isum+z_i{T}(lidx,:)];
      z2dum=[z_sumdum(1,:)*z_sumdum(1,:).' z_sumdum(2,:)*z_sumdum(2,:).'];
      N_c=[clust{clustidx,T}.numtemplates clust{nclustidx,T}.numtemplates]+[-1 1];
      ncost=((N_c.^2./(N_c-1))-z2dum./(fitlen*(N_c-1)));
      cVal=pcost-sum(ncost);
      
      %accept/reject cluster change
      if cVal>0
        reject=0;
      else
        if rand<exp(1/Tfac(T)*cVal)
          reject=0;
        else
          reject=1;
        end
      end
      
      if ~reject
        %update everything
        clustID(lidx,T)=nclustidx;
        clust{clustidx,T}.linidx=clust{clustidx,T}.linidx(~relidx);
        clust{nclustidx,T}.linidx=[clust{nclustidx,T}.linidx lidx];
        clust{clustidx,T}.z_isum=z_sumdum(1,:);
        clust{nclustidx,T}.z_isum=z_sumdum(2,:);
        clust{clustidx,T}.tempcost=ncost(1);
        clust{nclustidx,T}.tempcost=ncost(2);
        clust{clustidx,T}.numtemplates=clust{clustidx,T}.numtemplates-1;
        clust{nclustidx,T}.numtemplates=clust{nclustidx,T}.numtemplates+1;
        clust{nclustidx,T}.trl=[clust{nclustidx,T}.trl clust{clustidx,T}.trl(relidx)];
        clust{nclustidx,T}.tidx=[clust{nclustidx,T}.tidx clust{clustidx,T}.tidx(relidx)];
        clust{clustidx,T}.trl=clust{clustidx,T}.trl(~relidx);
        clust{clustidx,T}.tidx=clust{clustidx,T}.tidx(~relidx);        
        D(T)=D(T)-cVal;
        
        if mincost(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc{T};
            tclustID=clustID(:,T);
            tclust=clust(:,T);
            for nn=1:numclust
              clustnumel(nn)=tclust{nn}.numtemplates;
            end
          end
          mincost(T)=D(T);
          cc(T)=0;
        else
          cc(T)=cc(T)+1;
        end
      else
        rejcount(1,T)=rejcount(1,T)+1;
      end
      
        
     
    end
    

    
    
  end
  
  totcost(:,iter)=D;
  
  % Every couple of iterations, try a state swap between two random
  % temperatures
  if swapcount >= .5e3/(nPT-1);
    Tswap=randperm(nPT-1,1);
    Tswap=[Tswap Tswap+1];
    pswap=exp(1/konstant*diff(1./Tfac(Tswap))*diff(D(Tswap)));
    if rand < pswap;
      D(fliplr(Tswap))=D(Tswap);
      comcost(fliplr(Tswap))=comcost((Tswap));
      icomcost(:,fliplr(Tswap))=icomcost(:,(Tswap));
      loc(fliplr(Tswap))=loc(Tswap);
      z_i(fliplr(Tswap))=z_i(Tswap);
      clust(:,fliplr(Tswap))=clust(:,(Tswap));
      clustID(:,fliplr(Tswap))=clustID(:,(Tswap));
      
    end
    swapcount=0;
  end
  
  
  DTfac=[D; Tfac; cc];
  msg=sprintf([' # iterations: %d/%d\n cc =\n cost =           Tfac =         cc =\n' repmat('  %E     %1.4E  %8d\n',1, nPT) '\n Best clustering:\n ' repmat('%6d  ',1,numclust) '\n'], [iter numIt DTfac(:)' clustnumel]);
  if verbose
    fprintf([reverseStr, msg]);
  end
  reverseStr = repmat(sprintf('\b'), 1, length(msg));

end




totcost=totcost(:,1:iter);

% append the cost to previous run if possible
try
  cfg.totcost=[cfg.totcost; totcost.'];
catch
  cfg.totcost=totcost.';
end
cfg.finalcost=totcost(:,end);
cfg.cc=cc;
cfg.Tfac=(Tfac);
cfg.comcost=comcost;
cfg.icomcost=icomcost;
cfg.cm=cm;
try
cfg.bestloc=tloc;
catch
  cfg.bestloc=nan;
  tloc=loc{end};
end
try
cfg.bestclustID=tclustID;
catch
  cfg.bestclustID=clustID(:,end);
end

try
cfg.bestclust=tclust;
catch
  cfg.bestclust=clust(:,end);
end

cfg.clust=clust;
cfg.mincost=mincostTot;
cfg.loc=loc;

if ~exist('tclust','var')
  if numclust>1
  warning('no improvement in clustering found')
  end
  tclust=clust(:,nPT);
end

cfg.best_s=nan(fitlen,numclust);
cfg.best_z=nan(fitlen,numclust);
cfg.costdistr=cell(1,numclust);
for n=1:numclust
  
  dum_s=nan(tclust{n}.numtemplates,fitlen);
  for k=1:tclust{n}.numtemplates
    trl=tclust{n}.trl(k);
    tidx=tclust{n}.tidx(k);
    dum_s(k,:)=dat(trl,tloc(trl,tidx):tloc(trl,tidx)+fitlen-1);
  end
  cfg.best_s(:,n)=mean(dum_s,1);
  cfg.best_z(:,n)=mean(bsxfun(@rdivide,bsxfun(@minus,dum_s,mean(dum_s,2)),std(dum_s,1,2)),1);
  
  cfg.costdistr{n}=cost_i(dum_s);
end
cfg=orderfields(cfg);






% cleanup
fieldNamesdum=fieldnames(cfg);
[~, fieldNameIdx]=sort(lower(fieldNamesdum));
cfg=orderfields(cfg,fieldNameIdx);

if ~fullOutput
  cleanFields={'cc','clust', 'loc', 'finalcost', 'costdistr', 'comcost', 'cm','icomcost'};
  for nn=1:numel(cleanFields)
  try
    cfg=rmfield(cfg,cleanFields{nn});
  end
  end
  
  % undersample totcost
  if max(size(cfg.totcost))>2e3
    cfg.totcost_end=cfg.totcost(end-2e3:end,:);
    cfg.totcost_undSamp=cfg.totcost(1:1e2:end,:);
    cfg=rmfield(cfg,'totcost');
  end
end

return




function p = randval(n)
if numel(n)<2
  n=[n 1];
end
[~,p] = max(rand(n(1),n(2)));



function [loc,numWin]=initloc(guard,fitlen,data,numWin)
len=size(data);
dum=1:2*guard:max((len(2)-fitlen-guard),1);

if nargin<4
  numWin=numel(dum);
else
  if numWin>numel(dum)
    numWin=numel(dum);
  end
end



loc=nan(size(data,1),numWin+1);
for n=1:len(1)
  dum2=dum+randval([guard, numel(dum)])-1;
  dum2=min(dum2,len(2)-fitlen+1);
  sel=randperm(numel(dum2),numWin);
  loc(n,:)=[dum2(sel) size(data,2)+guard-fitlen];
end


function D_i=cost_i(z_s)

mu=mean(z_s,1);
z_s=bsxfun(@minus,z_s,mu).^2;
D_i=(z_s);

function [d cm d1]=com(y,ratio)
N=numel(y);
x=1:N;
yp=y-min(y);
cm=yp*x'/sum(yp);
d1=(cm-(x(end)+x(1))/2);
d=ratio*sumsqr(d1)/N;

function s2=sumsqr(x)
  selFinite = isfinite(x);
  x = x(selFinite);
  x2 = x.*x;
  s2 = sum(x2(:));
  
  
function sx=nanstd(x,varargin)
sel=~isnan(x);
sx=std(x(sel),varargin{1});