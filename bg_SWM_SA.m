function [cfg]=bg_SWM_SA(cfg, dat)
% [cfg]=bg_SWM(cfg, dat)
%
% Sliding Window Matching algorithm for detecting consistent, reocurring
% shapes in a data set. Using simulated annealing. Useful when minumum is
% found by bg_SWM, but temperature was still to high, and therefore
% convergence not yet complete.
%
% %%%%%%%
% INPUT:
% %%%%%%%
% cfg:  a structure that contains all parameters
% dat:  optional the data on which you want to apply the algorithm
%       (not neccesary if cfg.fname and cfg.varname are present)
%       Dimensions: M x N
%       M: trials
%       N: time points
%
% Fields in cfg that are required:
% .winLen:    the length of the sliding windows in timestamp units.
%             Note: This greatly affects the size of the found shape. I.e.
%             larger winLens tend to converge to shapes with lower-frequency
%             contents
% OR
% .winLenFac: Determines .winLen based on the peak in the 1/f-corrected
%             power spectrum of the data. If .winLen is omitted it will be
%             set to .winLenFac*T. Where T is the period corresponding to
%             the peak in the power spectrum.
%
% Optional fields (they all have a default value):
% CORE PARAMETERS
% .fname:     path and filename of the .mat-file that contains the data on
%             which the algorithms must be applied.
%             Note: if "dat" is given as input, this is set to 'function
%             input'
% .varname:   the name of the variable within .fname tha actually contains
%             the data.
%             Note: if "dat" is given as input, this is set to 'N/A'
% .numIt:     the number of iterations the algorithm will run through.
%             (default = 1e4)
% .guard:     the minimal space between the starting positions of two
%             succesive windows. Also determines the number of windows.
%             (default = .winLen/1.5)
% .numWin:    number of windows to fit per trial. Limited by length of
%             data, .winLen and .guard;
%             (default = max)
% .Tfac:      Temperature at which the SA algorithm starts.
%             (default = 1e-4)
%
% Note: good values for .Tfac depend on the scaling of your
% data and the signal to noise level
%
% MISC PARAMETERS
% .numClust:  number of shapes (=clusters) to find
%             (default = 1)
% .fs:        Sampling rate of dataset. Needed when using .Fbp or .Fbs
%             below
% .Fbp:       Bandpass frequency range. Usage: [hp lp] (rowvector; multiple
%             rows will apply multiple filters)
%             When .Fbp is supplied, the function first applies a bandpass
%             filter to the data before performing the SWM
% .Fbs:       Bandstop frequency range. Usage: [lp hp] (rowvector; multiple
%             rows will apply multiple filters)
%             When .Fbs is supplied, the function first applies a bandstop
%             filter to the data before performing the SWM
% .Fhp:       High-Pass cut-off frequency. use as above.
% .FhpFac:    High-Pass, but similar to .winLenFac, this applies filtering
%             relative to the peak frequency in the data. E.g. when the
%             main component is 10 Hz and .FhpFac=.5 it will apply 5 Hz HP
%             filtering.
% .Flp:       Low-Pass cut-off frequency. use as above.
%
% FLAGS
% .fullOutput:Flag to determine wheter the function should output the full
%             trajectory of the cost function, or just a undersampled
%             version of it. Useful for debugging or checking the effect of
%             chaning the temperature parameters
%             (default = 0)
% .verbose:   Determines whether the progression of the cost functions for
%             every temperature is printed to the screen. Useful if the
%             function is called directly. Not useful if the function is
%             sent as a batch job to a computer cluster that saves all
%             printed output. (Unneccesary large files are the result).
%             (default = 1)
% .dispPlot:  Determine wheter the cost function trajectory as well as
%             the mean shape is plotted.
%             i.e. a more visual realtive of .verbose.
%             (default = equal to .verbose)
%
% %%%%%%%
% OUTPUT:
% %%%%%%%
% cfg:	a structure containing the used parameters (i.e. the supplied and
% default values used during the algorithm run)
%
% Important fields:
% .best_s:        Contains the best mean shape for every cluster.
% .best_z:        Contains the best mean shape for every cluster; but after
%                 individual z-scoring (the algorithm uses this for calculating
%                 the cost function.
% .best_clust:    Contains the clustering state with lowest Cost.
% .best_clustID:  Contains the cluster IDs for every window for
%                 convenience. (No need to dig into best_clust).
% .costTotal:     The trajectories of the cost function for every temperature
%                 seperately. Only given when cfg.fullOutput=1.
% .totcost_unSamp:As above, but undersampled by a factor 100 to save
%                 diskspace. given when cfg.fullOutput=0;
% .costTotal_end: The cost values, but only for the final 2000 iterations.
%                 Given when cfg.fullOutput=0;
%
%
% Bart Gips 2014

%% check validity of input fields
validInp={'best_s';'best_z';'best_clust';'best_clustID';'best_loc';'clust';...
  'cm';'costCoM';'costDistr';'dispPlot';'Fbs';'Fbp';'Fhp';'FhpFac';'Flp';'costFinal';...
  'fname';'fs';'fullOutput';'guard';'costCoM_i';'loc';'costMin';...
  'numSelTemplates';'numClust';'numIt';'numIter';'numTemplates';'numWin';'ratio';'Tfac';...
  'costTotal';'costTotal_end';'costTotal_undSamp';'varname';'verbose';'winLen';'winLenFac';};
inpFields=fieldnames(cfg);

if any(~ismember(inpFields,validInp))
  badInp=inpFields(~ismember(inpFields,validInp));
  warning(sprintf(['Some fields in input cfg structure are invalid and are ignored:\n' repmat('  .%s\n',1,numel(badInp))],badInp{:}));
end



%% loading data
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

nanSel=isnan(dat);
nanFlag=any(nanSel(:));
if nanFlag
  warning(['Data contains ' num2str(round(sum(nanSel(:))/numel(nanSel)*100)) '% NaNs. Correct convergence is not guaranteed.'])
end

%% determine winLen
if isfield(cfg,'winLen')
  winLen=cfg.winLen;
else
  warning('winLen is not defined. Determining optimal winLen from 1/f corrected power spectrum')
  dat(nanSel)=0;
  nfft=2^nextpow2(size(dat,2))*4;
  spec=fft(diff(dat,1,2)',nfft);
  spec=spec(1:nfft/2+1,:);
  [~,midx]=max(mean(abs(spec).^2,2));
  shapeLen=nfft/midx;
  if size(dat,2)>200*shapeLen %multitaper to reduce noise
    clear spec
    mtLen=2^nextpow2(50*shapeLen)+1;
    taper=hanning(mtLen);
    reshapedum=ceil(numel(dat)/mtLen);
    datdum=zeros(size(dat,1),reshapedum*mtLen);
    datdum(:,1:size(dat,2))=dat;
    datdum=reshape(datdum.',mtLen,reshapedum).';
    datdum=bsxfun(@times,datdum,taper.');
    nfft=2^nextpow2(mtLen);
    spec=fft(diff(datdum,1,2)',nfft);
    spec=spec(1:nfft/2+1,:);
    [~,midx]=max(mean(abs(spec).^2,2));
    shapeLen=nfft/midx;
  end
  if isfield(cfg,'winLenFac')
    winLen=round(cfg.winLenFac*shapeLen);
  else
    winLen=round(1.5*shapeLen);
  end
  cfg.winLen=winLen;
end


%% preprocessing (filtering)
if any(isfield(cfg,{'Fbp','Fbs','Fhp','Flp','FhpFac'}))
  disp('Applying frequency filters...')
  filttype='fir';
  filttype='but';
  % first temporarily change NaNs to zero's such that filtering works
  if nanFlag
    dat(nanSel)=0;
    warning('Temporarily changing NaNs to zeros to make filtering possible');
  end
end



if isfield(cfg,'Fbp')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. Bandpass filter is not possible without cfg.fs')
  end
  Fbp=cfg.Fbp;
  if iscolumn(Fbp)
    Fbp=Fbp';
  end
  for band=1:size(Fbp,1)
    dat=ft_preproc_bandpassfilter(dat, fs, Fbp(band,:),[],filttype);
  end
end

if isfield(cfg,'Fbs')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. Bandstop filter is not possible without cfg.fs')
  end
  Fbs=cfg.Fbs;
  if iscolumn(Fbs)
    Fbs=Fbs';
  end
  for band=1:size(Fbs,1)
    dat=ft_preproc_bandstopfilter(dat, fs, Fbs(band,:),[],filttype);
  end
end

if isfield(cfg,'Fhp')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. High-pass filter is not possible without cfg.fs')
  end
  Fhp=cfg.Fhp;
  Fhp=Fhp(:);
  for freq=1:size(Fhp,1)
    dat=ft_preproc_highpassfilter(dat, fs, Fhp(freq),[],filttype);
  end
end

if isfield(cfg,'FhpFac')
  if exist('Fhp','var')
    warning('Applying double HP filtering. Based on both .Fhp and .FhpFac. This may not be desired.')
  end
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. High-pass filter is not possible without cfg.fs')
  end
  if ~exist('shapeLen','var')
    dat(nanSel)=0;
    nfft=2^nextpow2(size(dat,2))*4;
    spec=fft(diff(dat,1,2)',nfft);
    spec=spec(1:nfft/2+1,:);
    [~,midx]=max(mean(abs(spec).^2,2));
    shapeLen=nfft/midx;
    if size(dat,2)>200*shapeLen %multitaper to reduce noise
      clear spec
      mtLen=round(100*shapeLen);
      taper=hanning(mtLen);
      reshapedum=ceil(numel(dat)/mtLen);
      datdum=zeros(size(dat,1),reshapedum*mtLen);
      datdum(:,1:size(dat,2))=dat;
      datdum=reshape(datdum.',mtLen,reshapedum).';
      datdum=bsxfun(@times,datdum,taper.');
      nfft=2^nextpow2(mtLen)*4;
      spec=fft(diff(datdum,1,2)',nfft);
      spec=spec(1:nfft/2+1,:);
      [~,midx]=max(mean(abs(spec).^2,2));
      shapeLen=nfft/midx;
    end
  end
  Fhp=1/(shapeLen/fs)*cfg.FhpFac;
  Fhp=Fhp(:);
  for freq=1:size(Fhp,1)
    dat=ft_preproc_highpassfilter(dat, fs, Fhp(freq),[],filttype);
  end
end

if isfield(cfg,'Flp')
  if isfield(cfg,'fs')
    fs=cfg.fs;
  else
    error('Sampling rate missing. Low-pass filter is not possible without cfg.fs')
  end
  Flp=cfg.Flp;
  Flp=Flp(:);
  for freq=1:size(Flp,1)
    dat=ft_preproc_lowpassfilter(dat, fs, Flp(freq),[],filttype);
  end
end



%%

% change missing data back to NaNs
dat(nanSel)=nan;

if isfield(cfg,'guard')
  guard=cfg.guard;
else
  guard=round(winLen/1.5);
  cfg.guard=guard;
end


if isfield(cfg,'numIt')
  numIt=cfg.numIt;
else
  numIt=1e4;
  cfg.numIt=numIt;
end

% Do not use CoM weighting
if isfield(cfg,'ratio')
  %   ratio=cfg.ratio;
  ratio=0;
  cfg.ratio=ratio;
else
  ratio=0;
  cfg.ratio=ratio;
end


if isfield(cfg,'Tfac') %starting temperature
  Tfac=min(cfg.Tfac); %only take one temperature
else
  Tfac=10^-4;
  cfg.Tfac=Tfac;
end

nPT=numel(Tfac);
cfg.nPT=nPT;

if isfield(cfg,'numWin')
  numWin=cfg.numWin;
else
  numWin=0;
end

% initialize/read window locations
if isfield(cfg,'best_loc')
  loc=cfg.best_loc;
elseif isfield(cfg,'loc')
  loc=cfg.loc;
  if iscell(loc);
    loc=loc{1};
  end
  if size(loc,1)~=size(dat,1)
    error('location is not correct')
  end
  
else
  warning('location of windows missing (cfg.loc or cfg.best_loc); generating random initial window placement')
  if numWin
    [loc, numWin]=initloc(guard,winLen,dat,numWin);
  else
    [loc, numWin]=initloc(guard,winLen,dat);
  end
end


if numWin
  loc=loc(:,1:numWin);
else
  if var(loc(:,end))==0 %final column of loc is only a guard, not a real window
    loc=loc(:,1:end-1);
  end
  numWin=size(loc,2);
  cfg.numWin=numWin;
end



if isfield(cfg,'locSel')
  locSel=logical(cfg.locSel);
else
  locSel=true(size(loc));
end

if isfield(cfg,'maxShift')
  maxShift=cfg.maxShift;
else
  maxShift=round(.2*winLen);
  cfg.maxShift=maxShift;
end


if isfield(cfg,'numClust')
  numClust=cfg.numClust;
else
  cfg.numClust=1;
  numClust=1;
end

if isfield(cfg,'best_clust')
  clustInit=0;
  clust=cfg.best_clust;
  if isfield(cfg,'best_clustID')
    clustID=cfg.best_clustID;
  else
    clustID=nan(numTemplates,1);
    for n=1:numClust
      clustID(clust{n}.linIdx,1)=n;
    end
  end
else
  if isfield(cfg,'clust')
    clustInit=0;
    clust=cfg.clust;
  else
    clustInit=1;
  end
end

if isfield(cfg,'verbose')
  verbose=cfg.verbose;
else
  verbose=true;
end

if isfield(cfg,'dispPlot')
  dispPlot=cfg.dispPlot;
else
  dispPlot=verbose;
end


if isfield(cfg,'fullOutput')
  fullOutput=cfg.fullOutput;
else
  fullOutput=false;
  cfg.fullOutput=fullOutput;
end

%% Initialization

% find the number of sample vectors/templates
numTrl=size(dat,1);
numTemplates=numWin*numTrl;
cfg.numTemplates=numTemplates;
locSelInd=find(locSel);
numSelTemplates=numel(locSelInd);
cfg.numSelTemplates=numSelTemplates;
costCoM_i=nan(numTemplates,nPT); %CoM cost of individual windows

% construct sample vectors
z_i=nan(numTemplates,winLen);

if verbose
  fprintf('\nInitializing sliding windows and evaluating cost function...')
end
reverseStr='';



for n=1:numSelTemplates
  linidx=locSelInd(n);
  [trldum idxdum]=ind2sub([numTrl,numWin],linidx);
  s_i=dat(trldum,loc(trldum, idxdum):loc(trldum, idxdum)+winLen-1);
  z_i(linidx,:)=(s_i(:)-nanmean(s_i(:)))/nanstd(s_i(:),1);
  
  if ratio>0
    [costCoM_i(n,1) cm(n,1)]=com(z_i(linidx,:),ratio);
  else
    costCoM_i(n,1)=0;
    cm(n,1)=nan;
  end
  
  
  msg=sprintf('\n Template %d/%d', [n numSelTemplates]);
  if verbose
    fprintf([reverseStr, msg]);
  end
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
costCoM=nansum(costCoM_i(:,1));


if clustInit
  clust=cell(1);
  clustID=nan(numTemplates,1);
  clustID(:,1)=1;
  [trldum tidxdum]=ind2sub([numTrl, numWin],locSelInd);
  clust{1}.linIdx=locSelInd;
  clust{1}.trl=trldum;
  clust{1}.tidx=tidxdum;
  clust{1}.numTemplates=numel(locSelInd);
  
else %only reconstruct clust ID and initial cost
  if ~exist('clustID','var')
    clustID=nan(numTemplates,1);
    for n=1:numClust;
      clustID(clust{n}.linIdx,1)=n;
    end
  end
  
  %remove non-selected windows
  for n=1:numClust
    LIA=ismember(clust{n}.linIdx,locSelInd);
    clust{n}.linIdx=clust{n}.linIdx(LIA);
    clust{n}.trl=clust{n}.trl(LIA);
    clust{n}.tidx=clust{n}.tidx(LIA);
    clust{n}.numTemplates=sum(LIA);
  end
end



% initinalize the cost function
if verbose
  fprintf('\nInitializing cost matrices...')
end
D=zeros(1);
N=winLen;

for n=1:numClust
  N_c=clust{n}.numTemplates;
  z_isum=nansum(z_i(clust{n}.linIdx,:),1);
  clust{n}.z_isum=z_isum;
  if nanFlag
    % correct for NaNs
    clust{n}.noNanCount=sum(~isnan(z_i(clust{n}.linIdx,:)));
    nanFac=N_c./clust{n}.noNanCount;
    clust{n}.tempCost=(N_c^2/(N_c-1))-((nanFac.^2.*z_isum)*z_isum.')/(N*(N_c-1));
  else
    clust{n}.tempCost=(N_c^2/(N_c-1))-(z_isum*z_isum.')/(N*(N_c-1));
  end
  
  D=D+clust{n}.tempCost;
end


if verbose
  fprintf('\nDone!\n')
end



% draw plot window if requested
if dispPlot
  colorOrder=zeros(nPT,3);
  colorOrder(:,1)=linspace(0,1,nPT);
  colorOrder(:,3)=linspace(1,0,nPT); %make low T blue, high T Red
  plotSelT=1;
  colorOrderShape=zeros(numel(plotSelT),3);
  colorOrderShape(:,1)=linspace(0,1,numel(plotSelT));
  colorOrderShape(:,3)=linspace(1,0,numel(plotSelT)); %make low T blue, high T Red
  
  hfig=figure;
  set(hfig,'position',[100 100 1000 600])
  set(gcf,'DefaultAxesColorOrder',colorOrder)
  subplot(1,2,1)
  set(gca,'NextPlot','replacechildren')
  xlabel('iteration')
  ylabel('Cost')
  plotLegend=1;
  subplot(1,2,2)
end

TchangeIterCount=0;
TchangeToken=0;
stepSz=floor(guard/2);
TchangeCheck=1e3;
iterPlot=1;
rejcount=zeros(2,nPT);
clustnumel=nan(1,numClust);
stopToken=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Main optimization loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter=1;
costMin=D;
mincostTot=min(D);
costTotal=nan(nPT,numIt);
costTotal(:,1)=costMin;

cc=zeros(1,nPT);

reverseStr='';
if verbose
  fprintf('\nFinding templates...\n')
end

while iter<numIt &&  ~stopToken
  iter=iter+1;
  iterPlot=iterPlot+1;
  TchangeIterCount=TchangeIterCount+1;
  rejcount(2,:)=rejcount(2,:)+1;
  
  %determine to do cluster switch or window shift of sample
  pshift=.5;
  if numClust==1 %only one cluster means: never change clusters
    pshift=2;
  end
  
  for T=1:nPT
    lidx=randval(numTemplates);
    
    clustidx=clustID(lidx,T);
    
    shift=rand<pshift;
    
    if shift %shift a window (no cluster swap)
      [trl tidx]=ind2sub([numTrl,numWin],lidx);
      
      pLoc=loc(trl,tidx);
      %       pLoc=loc{T}(trl,tidx);
      locChange=0;
      loopcount=0;
      while ~locChange && loopcount<11 %find a possible step; Brute force for maximally 10 times
        dir=((rand>0.5)-0.5)*2*randval(stepSz);
        nLoc=pLoc+dir;
        locChange= nLoc>0 && size(dat,2)-nLoc>=winLen;
        
        if locChange
          %check guard:
          otherWinSel=true(numWin,1);
          otherWinSel(tidx)=false;
          minDist=min(abs(loc(trl,otherWinSel)-nLoc));
          locChange= minDist >= guard;
        end
        
        loopcount=loopcount+1;
      end
      
      
      
      if locChange %valid step is found, now evaluate change in cost function
        
        N_c=clust{clustidx,T}.numTemplates;
        %       pcost=D(T);
        pcost=clust{clustidx,T}.tempCost;
        pcomcost=costCoM_i(lidx,T);
        
        s_dum=dat(trl,nLoc:nLoc+winLen-1);
        z_dum=(s_dum-nanmean(s_dum))/nanstd(s_dum,1);
        
        if nanFlag
          pZ=z_i(lidx,:);
          pZ(isnan(pZ))=0;
          nZ=z_dum;
          nZ(isnan(nZ))=0;
          z_sumdum=clust{clustidx,T}.z_isum-pZ+nZ;
          noNanCountdum=clust{clustidx,T}.noNanCount-isnan(z_i(lidx,:))+isnan(z_dum);
          nanFac=N_c./noNanCountdum;
          ncost=(N_c^2/(N_c-1))-((nanFac.^2.*z_sumdum)*z_sumdum.')/(winLen*(N_c-1));
        else
          z_sumdum=clust{clustidx,T}.z_isum-z_i(lidx,:)+z_dum;
          ncost=(N_c^2/(N_c-1))-(z_sumdum*z_sumdum.')/(winLen*(N_c-1));
        end
        
        
        
        
        if ratio>0
          [ncomcost ncm]=com(z_dum,ratio);
        else
          ncomcost=0;
          ncm=nan;
        end
        % change in cost function
        cVal=(pcost-ncost+pcomcost-ncomcost);
      else %could not find a valid window shift -> "reject", i.e. keep everything as-is
        cVal=-inf;
      end
      
      % Determine wheter to accept or reject window shift
      if cVal>0
        accept=true;
      else
        if rand<exp(1/Tfac(T)*cVal)
          accept=true;
        else
          accept=false;
        end
      end
      
      if accept
        %update everything with new values
        loc(trl,tidx)=nLoc;
        %         loc{T}(trl,tidx)=nLoc;
        clust{clustidx,T}.z_isum=z_sumdum;
        if nanFlag
          clust{clustidx,T}.noNanCount=noNanCountdum;
        end
        clust{clustidx,T}.tempCost=ncost;
        D(T)=D(T)-cVal;
        z_i(lidx,:)=z_dum;
        costCoM_i(lidx,T)=ncomcost;
        costCoM(T)=costCoM(T)+ncomcost-pcomcost;
        cm(lidx,T)=ncm;
        
        if costMin(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc;
            %             tloc=loc{T};
          end
          costMin(T)=D(T);
          cc(T)=0;
        else
          cc(T)=cc(T)+1;
        end
      else
        rejcount(1,T)=rejcount(1,T)+1;
      end
      
      
    else %cluster swap instead of window shift
      relidx=clust{clustidx,T}.linIdx==lidx;
      nclustidx=randval(numClust);
      while nclustidx == clustidx
        nclustidx=randval(numClust);
      end
      pcost=clust{clustidx,T}.tempCost+clust{nclustidx,T}.tempCost;
      
      N_c=[clust{clustidx,T}.numTemplates clust{nclustidx,T}.numTemplates]+[-1 1];
      
      
      if nanFlag
        cZ=z_i(lidx,:);
        cZ(isnan(cZ))=0;
        z_sumdum=[clust{clustidx,T}.z_isum-cZ; clust{nclustidx,T}.z_isum+cZ];
        
        noNanCountdum=[clust{clustidx,T}.noNanCount-isnan(z_i(lidx,:)); clust{nclustidx,T}.noNanCount+isnan(z_i(lidx,:))];
        nanFac=bsxfun(@rdivide,N_c.',noNanCountdum);
        z2dum=[(nanFac(1,:).^2.*z_sumdum(1,:))*z_sumdum(1,:).' (nanFac(2,:).^2.*z_sumdum(2,:))*z_sumdum(2,:).'];
      else
        z_sumdum=[clust{clustidx,T}.z_isum-z_i(lidx,:); clust{nclustidx,T}.z_isum+z_i(lidx,:)];
        z2dum=[z_sumdum(1,:)*z_sumdum(1,:).' z_sumdum(2,:)*z_sumdum(2,:).'];
      end
      
      ncost=((N_c.^2./(N_c-1))-z2dum./(winLen*(N_c-1)));
      cVal=pcost-sum(ncost);
      
      %accept/reject cluster change
      if cVal>0
        accept=true;
      else
        if rand<exp(1/Tfac(T)*cVal)
          accept=true;
        else
          accept=false;
        end
      end
      
      if accept
        %update everything
        clustID(lidx,T)=nclustidx;
        clust{clustidx,T}.linIdx=clust{clustidx,T}.linIdx(~relidx);
        clust{nclustidx,T}.linIdx=[clust{nclustidx,T}.linIdx lidx];
        clust{clustidx,T}.z_isum=z_sumdum(1,:);
        clust{nclustidx,T}.z_isum=z_sumdum(2,:);
        clust{clustidx,T}.tempCost=ncost(1);
        clust{nclustidx,T}.tempCost=ncost(2);
        clust{clustidx,T}.numTemplates=clust{clustidx,T}.numTemplates-1;
        clust{nclustidx,T}.numTemplates=clust{nclustidx,T}.numTemplates+1;
        clust{nclustidx,T}.trl=[clust{nclustidx,T}.trl clust{clustidx,T}.trl(relidx)];
        clust{nclustidx,T}.tidx=[clust{nclustidx,T}.tidx clust{clustidx,T}.tidx(relidx)];
        clust{clustidx,T}.trl=clust{clustidx,T}.trl(~relidx);
        clust{clustidx,T}.tidx=clust{clustidx,T}.tidx(~relidx);
        D(T)=D(T)-cVal;
        
        if costMin(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc;
            %             tloc=loc{T};
            tclustID=clustID(:,T);
            tclust=clust(:,T);
            for nn=1:numClust
              clustnumel(nn)=tclust{nn}.numTemplates;
            end
          end
          costMin(T)=D(T);
          cc(T)=0;
        else
          cc(T)=cc(T)+1;
        end
      else
        rejcount(1,T)=rejcount(1,T)+1;
      end
    end
  end
  
  costTotal(:,iter)=D;
  
  % Every couple of iterations, try to lower the temperature (simulated
  % annealing)
  if TchangeIterCount >= TchangeCheck
    if TchangeToken > 5 % after five unsuccesful Temperature changes decrease stepsize
      if stepSz==2
        stopToken=1;
      end
      stepSz=ceil(stepSz/2+.5);
      TchangeToken=0;
    else
      costDum=costTotal(:,iter-TchangeCheck+1:iter);
      [P,s,mu]=polyfit(1:TchangeCheck,costDum,1);
      ste = sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
      if P(1)+2*ste(1) >0
        Tfac=Tfac/2;
        TchangeToken=1;
      else
        TchangeToken=0;
      end
    end
    TchangeIterCount=0;
  end
  
  DTfac=[D; Tfac; cc];
  msg=sprintf([' # iterations: %d/%d\n\n cost =           Tfac =    LastAccept =\n' repmat('  %E     %1.4E    %8d\n',1, nPT) '\n Best clustering:\n ' repmat('%6d  ',1,numClust) '\n'], [iter numIt DTfac(:)' clustnumel]);
  if verbose
    fprintf([reverseStr, msg]);
  end
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
  
  if dispPlot && iterPlot>100 % if requested; plot intermediate progress
    iterPlot=0;
    current_figure(hfig)
    subplot(1,2,1)
    plotselIter=max(1,iter-10e3):iter;
    plot(plotselIter,costTotal(:,plotselIter)','linewidth',2)
    xlim([plotselIter(1) plotselIter(1)+10e3-1])
    %     if plotLegend
    hleg=legend(num2str(Tfac(:),'%1.2e'));
    set(get(hleg,'title'),'string','Tfac')
    %       plotLegend=0;
    %     end
    subplot(1,2,2)
    
    for TT=1:numel(plotSelT)
      plot(nanmean(z_i),'color',colorOrderShape(TT,:),'linewidth',2)
      hold on
    end
    hold off
    title('mean shape (lowest temperatures)')
    hleg2=legend(num2str(Tfac(plotSelT)','%1.2e'));
    set(get(hleg2,'title'),'string','Tfac')
    drawnow
  end
  
  
end

costTotal=costTotal(:,1:iter);

%% create final output structure
% append the cost to previous run if possible
try
  cfg.costTotal=[cfg.costTotal; costTotal.'];
catch
  cfg.costTotal=costTotal.';
end
cfg.costFinal=costTotal(:,end);
cfg.costMin=mincostTot;
% if needed, CoM cost
if ratio>0
  cfg.costCoM=costCoM;
  cfg.costCoM_i=costCoM_i;
  cfg.cm=cm;
end

% make best_loc, only if improvement has been found
try
  cfg.best_loc=tloc;
catch
  cfg.best_loc=nan;
  tloc=loc;
end
try
  cfg.best_clustID=tclustID;
catch
  cfg.best_clustID=clustID(:,end);
end

% make best clustering, only if improvement has been found
try
  cfg.best_clust=tclust;
catch
  cfg.best_clust=clust(:,end);
end

% storing clusterings and locations for all temperatures, s.t. a second
% call to the function can resume where the first left off.
cfg.clust=clust;
cfg.loc=loc;

if ~exist('tclust','var')
  if numClust>1
    warning('No improvement in clustering found.')
  end
  tclust=clust(:,nPT);
end

% calculating the final shape (note this is after preprocessing, so this
% might yield different results than bg_swm_extract)
cfg.best_s=nan(winLen,numClust);
cfg.best_z=nan(winLen,numClust);
cfg.costDistr=cell(1,numClust);
for n=1:numClust
  
  dum_s=nan(tclust{n}.numTemplates,winLen);
  for k=1:tclust{n}.numTemplates
    trl=tclust{n}.trl(k);
    tidx=tclust{n}.tidx(k);
    dum_s(k,:)=dat(trl,tloc(trl,tidx):tloc(trl,tidx)+winLen-1);
  end
  cfg.best_s(:,n)=nanmean(dum_s,1);
  cfg.best_z(:,n)=nanmean(bsxfun(@rdivide,bsxfun(@minus,dum_s,mean(dum_s,2)),std(dum_s,1,2)),1);
  
  cfg.costDistr{n}=cost_i(z_score(dum_s));
end

cfg.Tfac=Tfac;

% cleanup
fieldNamesdum=fieldnames(cfg);
[~, fieldNameIdx]=sort(lower(fieldNamesdum));
cfg=orderfields(cfg,fieldNameIdx);

if ~fullOutput
  cleanFields={'cc','clust', 'loc', 'costFinal', 'costDistr', 'costCoM', 'cm','costCoM_i'};
  for nn=1:numel(cleanFields)
    try
      cfg=rmfield(cfg,cleanFields{nn});
    end
  end
  
  % undersample costTotal
  if max(size(cfg.costTotal))>2e3
    cfg.costTotal_end=cfg.costTotal(end-2e3:end,:);
    cfg.costTotal_undSamp=cfg.costTotal(1:1e2:end,:);
    cfg=rmfield(cfg,'costTotal');
  end
end

% sort in alphabetical order
cfg=orderfields(cfg);

return


%% subfunctions

function p = randval(n)
% alternative for randperm(n,1). but can do multiple at once (=n(2))
if numel(n)<2
  n=[n 1];
end
[~,p] = max(rand(n(1),n(2)));

function [loc,numWin]=initloc(guard,winLen,data,numWin)
% initialize window locations
len=size(data);
dum=1:2*guard:max((len(2)-winLen-guard),1);

if nargin<4
  numWin=numel(dum);
else
  if numWin>numel(dum)
    numWin=numel(dum);
  end
end

loc=nan(size(data,1),numWin);
for n=1:len(1)
  dum2=dum+randval([guard, numel(dum)])-1;
  dum2=min(dum2,len(2)-winLen+1);
  sel=randperm(numel(dum2),numWin);
  loc(n,:)=[dum2(sel)];
end

function D_i=cost_i(z_s)
% calculate Error (i.e. difference from mean)
mu=nanmean(z_s,1);
z_s=bsxfun(@minus,z_s,mu).^2;
D_i=(z_s);

function [d cm d1]=com(y,ratio)
% calculate centre of mass cost
N=numel(y);
x=1:N;
yp=y-min(y);
cm=yp*x'/sum(yp);
d1=(cm-(x(end)+x(1))/2);
d=ratio*sumsqr(d1)/N;

function s2=sumsqr(x)
% sum of squares
selFinite = isfinite(x);
x = x(selFinite);
x2 = x.*x;
s2 = sum(x2(:));

function sx=nanstd(x,varargin)
% replacement for nanstd
sel=~isnan(x);
sx=std(x(sel),varargin{1});

function current_figure(h)
set(0,'CurrentFigure',h)

function z=z_score(x)
x=bsxfun(@minus,x,nanmean(x));
z=bsxfun(@rdivide,x,sqrt(nanvar(x)));
