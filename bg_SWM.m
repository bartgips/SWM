function [cfg]=bg_SWM(cfg, dat)
% [cfg]=bg_SWM(cfg, dat)
%
% Sliding Window Matching algorithm for detecting consistent, reocurring
% shapes in a data set.
%
% %%%%%%%
% INPUT:
% %%%%%%%
% cfg:  a structure that contains all parameters
% dat:  optional the data on which you want to apply the algorithm
%       (not neccesary if cfg.fname and cfg.varname are present)
%       Dimensions: M x N x [Other dims]
%       M: trials
%       N: time points, dimension along which the windows slide
%       Other dims: other dimensions of the data, the windows do not slide
%       along this dimension, but they are taken into account for matching
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
% .nPT:       number of parallel temperatures
%             (default = 1)
% .Tfac:      Temperatures used in the PT sampling algorithm. This overrides
%             .nPT to numel(.Tfac).
%             (default = logspace(-3,1,.nPT))
% .konstant:  Bolzmann constant used in determining whether two
%             temperatures should be exchanged.
%             (default = 1e3)
% .mask:      A mask indicating "forbidden" parts of the data. I.e.
%             immovable barriers where the sliding windows cannot go. Mask
%             should be a logical the same size as the data. 0's mean no
%             mask. I.e. 1's will indicate "forbidden" data. If the data
%             contains NaN's, these will be automatically masked.
%             (note: can be sparse)
% .outputFile:A filename to where the output should be written. If this is
%             empty, the output will only be written to the workspace.
%             (Like most Matlab functions do).
% .zscore:    Calculate the cost function based on z-scored windows or not.
%             This is usefule when dealing with shapes of varying
%             amplitude.
%             (default = 0)
% .normalize: Normalize data before finding shape by means of z-scoring
%             (only really affects magnitude of values of cost function)
% 
%
% Note: good values for .Tfac and .konstant depend on the scaling of your
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
% .kernel:    A kernel that weights the cost function across the sliding
%             dimension. E.g. you are searching for an evoked response that
%             is more consistent at the beginning than at the end. Than a
%             linear kernel of linspace(1,0,winLen) may be useful. 
%             Note: the kernel can also extend over the other dimensions if
%             you whish to diferentially weight different parts of the
%             windows.
%             (default = [1], I.e. a unity kernel or no kernel)
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
%             initialization of the algorithm and is therefore useful during
%             debugging.
% .verbose:   Determines whether the progression of the cost functions for
%             every temperature is printed to the screen. Useful if the
%             function is called directly. Not useful if the function is
%             sent as a batch job to a computer cluster that saves all
%             printed output. (Unneccesary large files are the result).
%             (default = 1)
% .dispPlot:  Determine wheter the cost function trajectories as well as
%             the mean shapes with the lowest 3 temperatures is plotted.
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
% .costTotal_end:   The cost values, but only for the final 2000 iterations.
%                 Given when cfg.fullOutput=0;
%
%
% Bart Gips 2014

%% check validity of input fields
validInp={'best_s';'best_z';'best_clust';'best_clustID';'best_loc';'clust';...
  'costDistr';'costMin';'costFinal';'costTotal';...
  'costTotal_end';'costTotal_undSamp';'debug';'dispPlot';'Fbs';...
  'Fbp';'Fhp';'FhpFac';'Flp';'fname';'fs';'fullOutput';'guard';'kernel';...
  'konstant';'loc';'mask';'normalize';'nPT';'numClust';'numIt';'numIter';'numTemplates';...
  'numWin';'outputFile';'Tfac';'varname';'verbose';'winLen';'winLenFac';'zscore'};
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
  clear dum
  fInputFlag=false;
else
  cfg.fname='function input';
  cfg.varname='N/A';
  fInputFlag=true;
end

if iscolumn(dat);
  dat=dat.';
end

sz=size(dat);
if numel(sz)<3
  sz(3)=1;
  multiDim=0;
else
  multiDim=1;
end

nanSel=(isnan(dat));
nanFlag=any(nanSel(:));
if nanFlag
  warning(['Data contains ' num2str(round(sum(nanSel(:))/numel(nanSel)*100)) '% NaNs. Correct convergence is not guaranteed. Applying Mask..'])
  if isfield(cfg,'mask')
    mask=cfg.mask;
    mask=logical(mask+nanSel);
  else
    mask=nanSel;
  end
  cfg.mask=mask;
end

if isfield(cfg,'outputFile')
  outputFile=cfg.outputFile;
  outputFlag=true;
  % try writing to file
  try
    save(outputFile,'cfg');
  catch
    error(sprintf(['Unable to write to file:\n ' outputFile]));
  end
else
  outputFlag=false;
end

if isfield(cfg,'zscore')
  zscoreFlag=cfg.zscore;
else
  zscoreFlag=0;
end

if isfield(cfg,'normalize')
  normalize=cfg.normalize;
  if normalize && ~zscoreFlag
    datMean=nanmean(dat(:));
    datStd=nanstd(dat(:),1);
    dat=(dat-datMean)/datStd;
  else
    datMean=0;
    datStd=1;
  end
else
  datMean=0;
  datStd=1;
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
    for n=1:prod(sz(3:end))
      dat(:,:,n)=ft_preproc_bandpassfilter(dat(:,:,n), fs, Fbp(band,:),[],filttype);
    end
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
    for n=1:prod(sz(3:end))
      dat(:,:,n)=ft_preproc_bandstopfilter(dat(:,:,n), fs, Fbs(band,:),[],filttype);
    end
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
    for n=1:prod(sz(3:end))
      dat(:,:,n)=ft_preproc_highpassfilter(dat(:,:,n), fs, Fhp(freq),[],filttype);
    end
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
  cfg.Fhp=Fhp;
  for freq=1:size(Fhp,1)
    for n=1:prod(sz(3:end))
      dat(:,:,n)=ft_preproc_highpassfilter(dat(:,:,n), fs, Fhp(freq),[],filttype);
    end
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
    for n=1:prod(sz(3:end))
      dat(:,:,n)=ft_preproc_lowpassfilter(dat(:,:,n), fs, Flp(freq),[],filttype);
    end
  end
end



%%
% dat : mxn : m trials, n time points
% 2*guard should be bigger than winLen


% change missing data back to NaNs
dat(nanSel)=nan;

if isfield(cfg,'guard')
  guard=cfg.guard;
else
  guard=round(winLen/2);
  cfg.guard=guard;
end


if isfield(cfg,'numIt')
  numIt=cfg.numIt;
else
  numIt=1e4;
  cfg.numIt=numIt;
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
  Tfac=sort(cfg.Tfac); %sort it. The temperature swapping step requires that neighbouring temperatures are closest
  nPT=numel(Tfac);
  cfg.nPT=nPT;
else
  if nPT>1
    Tfac=logspace(-4,0,nPT);
  else
    Tfac=1e-4;
  end
  cfg.Tfac=Tfac;
end

if isfield(cfg,'numWin')
  numWin=cfg.numWin;
else
  numWin=0;
end

% initialize/read window locations
if isfield(cfg,'loc')
  loc=cfg.loc;
  if numel(loc)~= nPT || size(loc{1},1)~=size(dat,1)
    error('location is not correct')
  end
  if isfield(cfg,'best_loc')
    loc{1}=cfg.best_loc;
  end
else
  loc=cell(1,nPT);
  locStart=1;
  if isfield(cfg,'best_loc')
    loc{1}=cfg.best_loc;
    locStart=2;
  end
  for n=locStart:nPT
    if ~debug || n<2
      if numWin
        [loc{n}, numWin]=initloc(guard,winLen,dat,numWin);
      else
        [loc{n}, numWin]=initloc(guard,winLen,dat);
      end
    else
      loc{n}=loc{1};
    end
  end
  tloc=loc{end};
  cfg.numWin=numWin;
end

% prune window locations if mask is present
if isfield(cfg,'mask')
  maskFlag=true;
  mask=cfg.mask;
  numWinNew=0;
  nanMask=zeros(size(mask));
  nanMask(logical(mask))=nan;
  dat=bsxfun(@plus,dat,nanMask);
  for T=1:nPT
    if T==1 || ~debug
      for n=1:size(mask,1)
        for k=1:numWin
          selDum=loc{T}(n,k):min(loc{T}(n,k)+min(winLen,guard)-1,size(dat,2));
          if ~isnan(loc{T}(n,k)) && any(mask(n,selDum))
            loc{T}(n,k)=nan; % remove locs by making them into NaNs (to keep matrix dimensions of loc)
          end
        end
      end
      loc{T}=sort(loc{T},2);
      if numWinNew<numWin
        numWinNewDum=find(mean(isnan(loc{T}))>.95,1)-1;
        if isempty(numWinNewDum)
          numWinNew=numWin;
        else
          numWinNew=max(numWinNewDum,numWinNew);
        end
      end
    else
      loc{T}=loc{1};
    end
  end
  
  if numWinNew<numWin % check whether loc matrix can be made smaller
    numWin=numWinNew;
    cfg.numWin=numWin;
    for T=1:nPT
      loc{T}=loc{T}(:,1:numWin);
    end
  end
  
  if numWin<1
    error('Mask inhibits placement of all windows')
  end
else
  maskFlag=false;
end

if isfield(cfg,'kernel');
  kernel=cfg.kernel;
  if isrow(kernel)
    kernel=kernel(:);
  end
  kernel=kernel/sum(kernel(:))*numel(kernel); %normalize kernel to have comparable cost values with or without kernel
else
  kernel=1;
end


if isfield(cfg,'numClust')
  numClust=cfg.numClust;
else
  cfg.numClust=1;
  numClust=1;
end

if isfield(cfg,'clust')
  clustInit=0;
  clust=cfg.clust;
else
  clustInit=1;
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

% construct sample vectors
z_i=cell(1,nPT);
z_i(:)={nan([numTemplates,winLen,sz(3:end)])};

if verbose
  fprintf('\nInitializing sliding windows and evaluating cost function...')
end
reverseStr='';

for T=1:nPT
  if T<2 || ~debug
    for n=1:numTemplates
      [trldum idxdum]=ind2sub([numTrl,numWin],n);
      if maskFlag && isnan(loc{T}(trldum,idxdum))
        z_i{T}(n,:)=nan;
        nanFlag=1;
      else
        s_i=zeros([winLen,sz(3:end)]);
        s_i(:)=dat(trldum,loc{T}(trldum, idxdum):loc{T}(trldum, idxdum)+winLen-1,:);
        if zscoreFlag
          z_i{T}(n,:)=reshape(bsxfun(@times,kernel,(s_i-nanmean(s_i(:))))./nanstd(s_i(:),1),[],1);
        else
          z_i{T}(n,:)=reshape(bsxfun(@times,kernel,s_i),[],1);
        end
      end
      msg=sprintf('\n Temperature %d/%d \n Template %d/%d', [T nPT n numTemplates]);
      if verbose
        fprintf([reverseStr, msg]);
      end
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
      
      
    end
  else
    z_i{T}=z_i{1};
  end
end


if clustInit
  clust=cell(numClust,nPT);
  clustID=nan(numTemplates,T);
  tempperclust=round(numTemplates/numClust);
  for T=1:nPT
    % initialize clusters
    locidx=randperm(numTemplates);
    
    for clustidx=1:numClust
      
      if clustidx<numClust
        locidxdum=locidx(1:tempperclust);
        locidx=locidx(tempperclust+1:end);
      else
        locidxdum=locidx;
      end
      clustID(locidxdum,T)=clustidx;
      [trldum tidxdum]=ind2sub([numTrl, numWin],locidxdum);
      clust{clustidx,T}.linIdx=locidxdum;
      clust{clustidx,T}.trl=trldum;
      clust{clustidx,T}.tidx=tidxdum;
      clust{clustidx,T}.numTemplates=numel(locidxdum);
      
    end
  end
  
else %only reconstruct clust ID and initial cost
  clustID=nan(numTemplates,nPT);
  for T=1:nPT
    for n=1:numClust;
      clustID(clust{n,T}.linIdx,T)=n;
      %     clust{n}.z_isum=sum(z_i(clust{n}.linIdx,:),1);
    end
  end
end

if isfield(cfg,'best_clust')
  clust(:,nPT)=cfg.best_clust;
  for n=1:numClust
    clustID(clust{n,nPT}.linIdx,nPT)=n;
  end
end

% initinalize the cost function
if verbose
  fprintf('\nInitializing cost matrices...')
end
D=zeros(1,nPT);
N=winLen*prod(sz(3:end));
for T=1:nPT
  D(T)=0;
  for n=1:numClust
    N_c=clust{n,T}.numTemplates;
    z_isum=nansum(z_i{T}(clust{n,T}.linIdx,:),1);
    if zscoreFlag
      clust{n,T}.varianceFac=N_c;
    else
      clust{n,T}.varianceFac=N_c*nanmean(nanmean(z_i{T}(clust{n,T}.linIdx,:).^2,2));
    end
    varianceFac=clust{n,T}.varianceFac;
    clust{n,T}.z_isum=z_isum;
    if nanFlag
      % correct for NaNs
      clust{n,T}.noNanCount=sum(~isnan(z_i{T}(clust{n,T}.linIdx,:)));
      nanFac=N_c./clust{n,T}.noNanCount;
      clust{n,T}.tempCost=(varianceFac*N_c/(N_c-1))-((nanFac.^2.*z_isum)*z_isum.')/(N*(N_c-1));
    else
      clust{n,T}.tempCost=(varianceFac*N_c/(N_c-1))-(z_isum*z_isum.')/(N*(N_c-1));
    end
    
    D(T)=D(T)+clust{n,T}.tempCost;
  end
end

if verbose
  fprintf('\nDone!\n')
end

% draw plot window if requested
if dispPlot
  colorOrder=zeros(nPT,3);
  colorOrder(:,1)=linspace(0,1,nPT);
  colorOrder(:,3)=linspace(1,0,nPT); %make low T blue, high T Red
  
  plotselT=1:min(3,nPT); % plot shapes of lowest 3 temperatures.
  colorOrderShape=zeros(numel(plotselT),3);
  colorOrderShape(:,1)=linspace(0,1,numel(plotselT));
  colorOrderShape(:,3)=linspace(1,0,numel(plotselT)); %make low T blue, high T Red
  
  hfig=figure;
  set(hfig,'position',[100 100 1000 600])  
  subplot(1,2,1)
  set(gca,'ColorOrder',flipud(colorOrder))
  set(gca,'NextPlot','replacechildren')
  xlabel('iteration')
  ylabel('Cost')
  plotLegend=1;
end

swapcount=0;
Tchangecount=0;
saveCount=0;
iterPlot=1;
rejcount=zeros(2,nPT);
clustnumel=nan(1,numClust);

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

while iter<numIt %&&  cc<cclim
  iter=iter+1;
  iterPlot=iterPlot+1;
  swapcount=swapcount+1;
  saveCount=saveCount+outputFlag;
  
  Tchangecount=Tchangecount+1;
  rejcount(2,:)=rejcount(2,:)+1;
  
  %determine to do cluster switch or window shift of sample
  pshift=.5;
  if numClust==1 %only one cluster means: never change clusters
    pshift=2;
  end
  
  for T=1:nPT
    lidx=randval(numTemplates);
    
    while maskFlag && isnan(loc{T}(lidx)) % do not try to move pruned locs
      lidx=randval(numTemplates);
    end
    clustidx=clustID(lidx,T);
    
    shift=rand<pshift;
    
    if shift %shift a window (no cluster swap)
      [trl tidx]=ind2sub([numTrl,numWin],lidx);
      
      pLoc=loc{T}(trl,tidx);
      locChange=0;
      loopcount=0;
      while ~locChange && loopcount<11 %find a possible step; Brute force for maximally 10 times
        dir=((rand>0.5)-0.5)*2*randval(floor(guard/2));
        nLoc=pLoc+dir;
        locChange= nLoc>0 && size(dat,2)-nLoc>=winLen;
        
        if locChange && numWin>1
          %check guard:
          otherWinSel=true(numWin,1);
          otherWinSel(tidx)=false;
          minDist=min(abs(loc{T}(trl,otherWinSel)-nLoc));
          locChange= minDist >= guard;
        end
        
        % also check for distance from mask
        if locChange && maskFlag
          selDum=nLoc:min(nLoc+min(winLen,guard)-1,size(dat,2));
          maskDum=sum(mask(trl,selDum))<1;
          locChange=maskDum;
        end
        loopcount=loopcount+1;
      end
      
      
      
      if locChange %valid step is found, now evaluate change in cost function
        
        N_c=clust{clustidx,T}.numTemplates;
        %       pcost=D(T);
        pcost=clust{clustidx,T}.tempCost;
        s_dum=zeros(winLen,sz(3:end));
        s_dum(:)=dat(trl,nLoc:nLoc+winLen-1,:);
        
        if zscoreFlag
          z_dum=bsxfun(@times,kernel,(s_dum-nanmean(s_dum(:)))/nanstd(s_dum(:),1));
          varianceFac=N_c;
        else
          z_dum=bsxfun(@times,kernel,s_dum);
          varianceFac=clust{clustidx,T}.varianceFac;
          pZ=z_i{T}(lidx,:);
          varianceFac=varianceFac-nanmean(pZ.^2)+nanmean(z_dum(:).^2);
        end
        
        
        
        if nanFlag
          pZ=z_i{T}(lidx,:);
          pZ(isnan(pZ))=0;
          nZ=z_dum;
          nZ(isnan(nZ))=0;
          z_sumdum=clust{clustidx,T}.z_isum-pZ+nZ(:)';
          noNanCountdum=clust{clustidx,T}.noNanCount-isnan(z_i{T}(lidx,:))+isnan(z_dum(:)');
          nanFac=N_c./noNanCountdum;
          ncost=(varianceFac*N_c/(N_c-1))-((nanFac.^2.*z_sumdum)*z_sumdum.')/(winLen*prod(sz(3:end))*(N_c-1));
        else
          z_sumdum=clust{clustidx,T}.z_isum-z_i{T}(lidx,:)+z_dum(:)';
          ncost=(varianceFac*N_c/(N_c-1))-(z_sumdum*z_sumdum.')/(winLen*prod(sz(3:end))*(N_c-1));
        end
        
        
        
        % change in cost function
        cVal=(pcost-ncost);
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
        cc(T)=0;
        %update everything with new values
        loc{T}(trl,tidx)=nLoc;
        clust{clustidx,T}.z_isum=z_sumdum;
        clust{clustidx,T}.varianceFac=varianceFac;
        if nanFlag
          clust{clustidx,T}.noNanCount=noNanCountdum;
        end
        clust{clustidx,T}.tempCost=ncost;
        D(T)=D(T)-cVal;
        z_i{T}(lidx,:)=z_dum(:);
        if costMin(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc{T};
          end
          costMin(T)=D(T);
          
        end
      else
        cc(T)=cc(T)+1;
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
      
      cZ=z_i{T}(lidx,:);
      if nanFlag
        
        cZ(isnan(cZ))=0;
        z_sumdum=[clust{clustidx,T}.z_isum-cZ; clust{nclustidx,T}.z_isum+cZ];
        
        noNanCountdum=[clust{clustidx,T}.noNanCount-isnan(z_i{T}(lidx,:)); clust{nclustidx,T}.noNanCount+isnan(z_i{T}(lidx,:))];
        nanFac=bsxfun(@rdivide,N_c.',noNanCountdum);
        z2dum=[(nanFac(1,:).^2.*z_sumdum(1,:))*z_sumdum(1,:).' (nanFac(2,:).^2.*z_sumdum(2,:))*z_sumdum(2,:).'];
      else
        z_sumdum=[clust{clustidx,T}.z_isum-cZ; clust{nclustidx,T}.z_isum+cZ];
        z2dum=[z_sumdum(1,:)*z_sumdum(1,:).' z_sumdum(2,:)*z_sumdum(2,:).'];
      end
      
      if zscoreFlag
        varianceFac=N_c;
      else
        varianceFac=[clust{clustidx,T}.varianceFac clust{nclustidx,T}.varianceFac];
        varianceFac=varianceFac+[-nanmean(cZ.^2) nanmean(cZ.^2)];
      end
      
      ncost=((varianceFac.*N_c./(N_c-1))-z2dum./(winLen*prod(sz(3:end))*(N_c-1)));
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
        cc(T)=0;
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
        clust{clustidx,T}.varianceFac=varianceFac(1);
        clust{nclustidx,T}.varianceFac=varianceFac(2);
        D(T)=D(T)-cVal;
        
        if costMin(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc{T};
            tclustID=clustID(:,T);
            tclust=clust(:,T);
            for nn=1:numClust
              clustnumel(nn)=tclust{nn}.numTemplates;
            end
          end
          costMin(T)=D(T);
        end
      else
        cc(T)=cc(T)+1;
        rejcount(1,T)=rejcount(1,T)+1;
      end
    end
  end
  
  costTotal(:,iter)=D;
  
  % Every couple of iterations, try a state swap between two randomly
  % chosen temperatures
  if nPT>1 && swapcount >= .5e3/(nPT);
    Tswap=randperm(nPT,2);
    %     Tswap=[Tswap Tswap+1];
    pswap=exp(1/konstant*diff(1./Tfac(Tswap))*diff(D(Tswap)));
    if rand < pswap;
      D(fliplr(Tswap))=D(Tswap);
      loc(fliplr(Tswap))=loc(Tswap);
      z_i(fliplr(Tswap))=z_i(Tswap);
      clust(:,fliplr(Tswap))=clust(:,(Tswap));
      clustID(:,fliplr(Tswap))=clustID(:,(Tswap));
    end
    swapcount=0;
  end
  
  DTfac=[D; Tfac; cc];
  if numClust>1
    msg=sprintf([' # iterations: %d/%d\n\n cost =           Tfac =    LastAccept =\n' repmat('  %E     %1.4E    %8d\n',1, nPT) '\n Best clustering:\n ' repmat('%6d  ',1,numClust) '\n'], [iter numIt DTfac(:)' clustnumel]);
  else
    msg=sprintf([' # iterations: %d/%d\n\n cost =           Tfac =    LastAccept =\n' repmat('  %E     %1.4E    %8d\n',1, nPT) '\n'], [iter numIt DTfac(:)']);
  end
  if verbose
    fprintf([reverseStr, msg]);
  end
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
  
  if dispPlot && iterPlot>50 % if requested; plot intermediate progress
    iterPlot=0;
    current_figure(hfig)
    subplot(1,2,1)
    plotselIter=max(1,iter-5e3):iter;
    h=plot(plotselIter,fliplr(costTotal(:,plotselIter)'),'linewidth',2);
    xlim([plotselIter(1) plotselIter(1)+5e3-1])
    if plotLegend
      if numel(Tfac)>5
        TfacSel=unique(round(linspace(1,numel(Tfac),5)));
        TfacDum=Tfac(:);
        hleg=legend(h(TfacSel),num2str(flipud(TfacDum(TfacSel)),'%1.2e'),'location','southwest');
      else
        hleg=legend(num2str(flipud(Tfac(:)),'%1.2e'),'location','southwest');
      end
      set(get(hleg,'title'),'string','Tfac')
      plotLegend=0;
    end
    
    if multiDim
      for n=1:min(numClust,3)
        subplot(min(numClust,3),2,sub2ind([2 min(numClust,3)],2,n))
        imDum=datMean+datStd*bsxfun(@rdivide,reshape(clust{n,1}.z_isum/clust{n,1}.numTemplates,winLen,[]),kernel)';
        imagesc(imDum,maxabs(imDum));
        h=colorbar;
        if zscoreFlag
          set(get(h,'ylabel'),'string','z-score')
        else
          set(get(h,'ylabel'),'string','mean Signal')
        end
        if numClust>1
          title(['mean shape (lowest temperature); cluster ' num2str(n) '/' num2str(numClust)])
        else
          title(['mean shape (lowest temperature)'])
        end
        ylabel('Concatenated other dimensions')
      end
    else
      if numClust<2
      for TT=1:numel(plotselT)
        plot(clust{n,datMean+datStd*plotselT(TT)}.z_isum/clust{n,plotselT(TT)}.numTemplates./kernel,'color',colorOrderShape(TT,:),'linewidth',2)
        hold on
      end
      hold off
      hleg2=legend(num2str(Tfac(plotselT)','%1.2e'));
      set(get(hleg2,'title'),'string','Tfac')
      title('mean shape (lowest temperatures)')
      if zscoreFlag
        ylabel('z-score')
      else
        ylabel('mean Signal')
      end
      else
        for TT=1:numel(plotselT)
          subplot(numel(plotselT),2,sub2ind([2 numel(plotselT)],2,TT))
          set(gca,'ColorOrder',lines(8))
          for n=1:min(numClust)
          plot(datMean+datStd*clust{n,plotselT(TT)}.z_isum/clust{n,plotselT(TT)}.numTemplates./kernel,'linewidth',2)
          hold all
          end
          hold off
          title(['mean shapes (Tfac =' num2str(Tfac(plotselT(TT)),'%1.2e') ')'])
        end
      end
    end
    xlabel('Sliding dimension')
    drawnow
  end
  
  if saveCount >1e3 || iter==numIt  % every 1000 iterations, save file to disk
    saveCount=0;
    %% create final output structure
    % append the cost to previous run if possible
    %     try
    %       cfg.costTotal=[cfg.costTotal; costTotal.'];
    %     catch
    cfg.costTotal=costTotal.';
    %     end
    cfg.costFinal=costTotal(:,end);
    cfg.costMin=mincostTot;
    
    % make best_loc, only if improvement has been found
    try
      cfg.best_loc=tloc;
    catch
      cfg.best_loc=nan;
      tloc=loc{1};
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
    cfg.best_s=nan([numClust,winLen,sz(3:end)]);
    cfg.best_z=cfg.best_s;
    cfg.costDistr=cell(1,numClust);
    for n=1:numClust
      
      dum_s=nan([tclust{n}.numTemplates,winLen, sz(3:end)]);
      for k=1:tclust{n}.numTemplates
        trl=tclust{n}.trl(k);
        tidx=tclust{n}.tidx(k);
        if ~isnan(tloc(trl,tidx))
          dum_s(k,:)=reshape(dat(trl,tloc(trl,tidx):tloc(trl,tidx)+winLen-1,:),[],1);
        end
      end
      cfg.best_s(n,:)=reshape(nanmean(dum_s,1),[],1);
      cfg.best_z(n,:)=nanmean(bsxfun(@rdivide,bsxfun(@minus,reshape(dum_s,size(dum_s,1),[]),mean(reshape(dum_s,size(dum_s,1),[]),2)),std(reshape(dum_s,size(dum_s,1),[]),1,2)),1);
      
      cfg.costDistr{n}=cost_i(z_score(dum_s));
    end
    % make clusters the last dimension
    cfg.best_s=datMean+datStd*permute(cfg.best_s,[2:ndims(cfg.best_s),1]);
    cfg.best_z=permute(cfg.best_z,[2:ndims(cfg.best_z),1]);
    % sort in alphabetical order
    cfg=orderfields(cfg);
    
    % cleanup
    fieldNamesdum=fieldnames(cfg);
    [~, fieldNameIdx]=sort(lower(fieldNamesdum));
    cfg=orderfields(cfg,fieldNameIdx);
    
    if ~fullOutput
      cleanFields={'cc','clust', 'costFinal', 'costDistr', 'costCoM', 'cm','costCoM_i'};
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
    
    if outputFlag
      save(outputFile,'cfg');
      if fInputFlag
        save(outputFile,'dat','-append');
      end
    end
    
  end
  
end

costTotal=costTotal(:,1:iter);



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
maxWin= floor(len(2)/guard)-1;
dum=guard/2:1.5*guard:max((len(2)-winLen),1);
if nargin<4
  numWin=maxWin;
else
  if numWin>maxWin
    numWin=maxWin;
  end
end

loc=nan(size(data,1),numWin);

nEmptySpace=len(2)-winLen-(numWin-1)*guard;
guards=[0:numWin-1]*guard;

for n=1:len(1)
  empty=sort(randperm(nEmptySpace,numWin));
  winStarts=empty+guards;
  
  loc(n,:)=sort(winStarts);
end

function D_i=cost_i(z_s)
% calculate Error (i.e. difference from mean)
mu=nanmean(z_s,1);
z_s=bsxfun(@minus,z_s,mu).^2;
D_i=(z_s);

function s2=sumsqr(x)
% sum of squares
selFinite = isfinite(x);
x = x(selFinite);
x2 = x.*x;
s2 = sum(x2(:));

function sx=nanstd(x,varargin)
% replacement for nanstd
sel=~isnan(x);
if nargin>1
  sx=std(x(sel),varargin{1});
else
  sx=std(x(sel));
end

function current_figure(h)
set(0,'CurrentFigure',h)

function z=z_score(x)
x=bsxfun(@minus,x,nanmean(x));
z=bsxfun(@rdivide,x,sqrt(nanvar(x)));
