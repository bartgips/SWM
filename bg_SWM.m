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
%       (not neccesary if cfg.fName and cfg.varName are present)
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
% .fName:     path and filename of the .mat-file that contains the data on
%             which the algorithms must be applied.
%             Note: if "dat" is given as input, this is set to 'function
%             input'
% .varName:   the name of the variable within .fName tha actually contains
%             the data.
%             Note: if "dat" is given as input, this is set to 'N/A'
% .numIt:     the number of iterations the algorithm will run through.
%             (default = 1e4)
% .guard:     the minimal space between the starting positions of two
%             succesive windows. Also determines the number of windows.
%             (default = .winLen/1.5)
% .winPerTrial:    number of windows to fit per trial. Limited by length of
%             data, .winLen and .guard;
%             (default = max)
% .nPT:       number of parallel temperatures
%             (default = 1)
% .Tfac:      Temperatures used in the PT sampling algorithm. This overrides
%             .nPT to numel(.Tfac).
%             (default = logspace(-3,1,.nPT))
% .konstant:  Bolzmann constant used in determining whether two
%             temperatures should be exchanged.
%             (default = 100, higher values make switches more likely)
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
%             (default = 1)
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
% .FoVInterval: The average number of times all windows are moved, before
%             the entire FoV is shifted. An FoV shift means shifting all
%             windows within a cluster at once. This makes sure the mean
%             shape is "centered" correctly inside the windows.
%             (default = 2; i.e. FoV-shift is attempted every 2*numWindows
%             iterations.)
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
% .costTotal_unSamp:As above, but undersampled by a factor 100 to save
%                 diskspace. given when cfg.fullOutput=0;
% .costTotal_end: The cost values, but only for the final 2000 iterations.
%                 Given when cfg.fullOutput=0;
%
%
% Bart Gips 2014

compatibilityFlag=verLessThan('matlab', '8.4');
%% check validity of input fields
validInp={'best_s';'best_z';'best_clust';'best_clustID';'best_loc';'clust';...
  'costDistr';'costMin';'costFinal';'costTotal';...
  'costTotal_end';'costTotal_undSamp';'debug';'dispPlot';'Fbs';...
  'Fbp';'Fhp';'FhpFac';'Flp';'fName';'fs';'FoVInterval';'fullOutput';'guard';'kernel';...
  'konstant';'loc';'mask';'normalize';'nPT';'numClust';'numIt';'numIter';'numWindows';...
  'winPerTrial';'outputFile';'stepSz';'Tfac';'varName';'verbose';'winLen';'winLenFac';'zscore'};

% fix case of inputs
try
  [~,cfg]=isfieldi(cfg,validInp);
catch
  error('cfg contains conflicting fields')
end

inpFields=fieldnames(cfg);

if any(~ismember(inpFields,validInp))
  badInp=inpFields(~ismember(inpFields,validInp));
  warning(sprintf(['Some fields in input cfg structure are invalid and are ignored:\n' repmat('  .%s\n',1,numel(badInp))],badInp{:}));
end

%% loading data
% load data from file instead of from function input
if nargin<2
  dum=load(cfg.fName, cfg.varName);
  eval(['dat=dum.' cfg.varName ';'])
  clear dum
  fInputFlag=false;
else
  cfg.fName='function input';
  cfg.varName='N/A';
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

if isfield(cfg,'verbose')
  verbose=cfg.verbose;
else
  verbose=true;
end


if isfield(cfg,'outputFile')
  outputFile=cfg.outputFile;
  outputFlag=true;
  [outputDir, outputFile]=fileparts(outputFile); % remove file extension if present
  outputFile=[fullfile(outputDir,outputFile) '.mat']; % adding .mat extension
  if verbose && exist([outputFile],'file')
    warning('outputFile already exists...')
    reply = input('Do you want to overwrite? Y/N [N]: ', 's');
    if isempty(reply)
      outputFlag = false;
    elseif strcmpi(reply,'n')
      outputFlag = false;
    elseif strcmpi(reply,'y')
      outputFlag = true;
    else
      warning('input not recognized, not overwriting')
      outputFlag = false;
    end
  end
  
  if outputFlag
    % try writing to file
    try
      save(outputFile,'cfg');
    catch
      error(sprintf(['Unable to write to file:\n ' outputFile]));
    end
  end
else
  outputFlag=false;
end

if isfield(cfg,'zscore')
  zscoreFlag=cfg.zscore;
else
  zscoreFlag=1;
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
  [~,mIdx]=max(mean(abs(spec).^2,2));
  shapeLen=nfft/mIdx;
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
    [~,mIdx]=max(mean(abs(spec).^2,2));
    shapeLen=nfft/mIdx;
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
    [~,mIdx]=max(mean(abs(spec).^2,2));
    shapeLen=nfft/mIdx;
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
      [~,mIdx]=max(mean(abs(spec).^2,2));
      shapeLen=nfft/mIdx;
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



%% Determining parameters and their defaults
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
  konstant=1e2;
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

if isfield(cfg,'winPerTrial')
  winPerTrial=cfg.winPerTrial;
else
  winPerTrial=0;
end

%read mask
if isfield(cfg,'mask')
  maskFlag=true;
  mask=cfg.mask;
else
  mask=[];
  maskFlag=false;
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
  if ~winPerTrial
    winPerTrial=size(loc{1},2);
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
      if winPerTrial
        [loc{n}, winPerTrial, stepSzdum]=initloc(guard,winLen,dat,winPerTrial,mask);
      else
        [loc{n}, winPerTrial, stepSzdum]=initloc(guard,winLen,dat,[],mask);
      end
    else
      loc{n}=loc{1};
    end
  end
  tloc=loc{end};
  cfg.winPerTrial=winPerTrial;
end

if isfield(cfg,'stepSz')
  stepSz=cfg.stepSz;
else
  if exist('stepSzdum','var') % stepsize determined by initloc
    stepSz=stepSzdum;
    clear stepSzdum;
  else
    stepSz=ceil((size(dat,2)-(winPerTrial-1)*guard-winLen)/winPerTrial);
  end
  cfg.stepSz=stepSz;
end

% prune window locations if mask is present
if maskFlag
  winPerTrialNew=0;
  nanMask=zeros(size(mask));
  nanMask(logical(mask))=nan;
  dat=bsxfun(@plus,dat,nanMask);
  for T=1:nPT
    if T==1 || ~debug
      for n=1:size(mask,1)
        for k=1:winPerTrial
          selDum=loc{T}(n,k):min(loc{T}(n,k)+min(winLen,guard)-1,size(dat,2));
          selDum=selDum(selDum>0 & selDum <= sz(2));
          if ~isnan(loc{T}(n,k)) && any(mask(n,selDum))
            loc{T}(n,k)=nan; % remove locs by making them into NaNs (to keep matrix dimensions of loc)
          end
        end
      end
      loc{T}=sort(loc{T},2);
      if winPerTrialNew<winPerTrial
        winPerTrialNewDum=find(mean(isnan(loc{T}))>.95,1)-1;
        if isempty(winPerTrialNewDum)
          winPerTrialNew=winPerTrial;
        else
          winPerTrialNew=max(winPerTrialNewDum,winPerTrialNew);
        end
      end
    else
      loc{T}=loc{1};
    end
  end
  
  if winPerTrialNew<winPerTrial % check whether loc matrix can be made smaller
    winPerTrial=winPerTrialNew;
    cfg.winPerTrial=winPerTrial;
    for T=1:nPT
      loc{T}=loc{T}(:,1:winPerTrial);
    end
  end
  
  if winPerTrial<1
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

validClustFields={'linIdx','trl','tIdx','numWindows','tempCost','z_isum',...
  'varianceFac'};
if isfield(cfg,'clust')
  clustInit=0;
  clust=cfg.clust;
  for n=1:numel(clust)
    [~,clust{n}]=isfieldi(clust{n},validClustFields);
  end
else
  clustInit=1;
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

if isfield(cfg,'FoVInterval')
  FoVInterval=cfg.FoVInterval;
else
  FoVInterval=1;
  cfg.FoVInterval=FoVInterval;
end

%% Initialization

% find the number of sample vectors/templates
numTrl=size(dat,1);
numWindows=winPerTrial*numTrl;
cfg.numWindows=numWindows;

% determining Tswap and FoV tunneling intervals based on numWindows and
% nPT.
% Shift all windows 5 times on average before attempting to swap Tfac
TswapInterval=5*numWindows/nPT;
% Shift all windows FoVInterval times on average before attempting to shift FoV
FoVShiftInterval=FoVInterval*numWindows;

% construct sample vectors
z_i=cell(1,nPT);
z_i(:)={nan([numWindows,winLen,sz(3:end)])};

if verbose
  fprintf('\nInitializing sliding windows and evaluating cost function...')
end
reverseStr='';

for T=1:nPT
  if T<2 || ~debug
      
    msg=sprintf('\n Temperature %d/%d', [T nPT]);
    if verbose
      fprintf([reverseStr, msg]);
    end
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    cfgExtr=[];
    cfgExtr.loc=loc{T};
    cfgExtr.winLen=cfg.winLen;
    cfgExtr.numWindows=cfg.numWindows;
    
    [s,z]=bg_swm_extract(cfgExtr,dat);
    
    if zscoreFlag
      z_i{T}=z;
    else
      z_i{T}=s;
    end
    
    if any(isnan(z_i{T}))
      nanFlag=1;
    end
  else
    z_i{T}=z_i{1};
  end
end


if clustInit
  clust=cell(numClust,nPT);
  clustID=nan(numWindows,T);
  tempperclust=round(numWindows/numClust);
  for T=1:nPT
    % initialize clusters
    locIdx=randperm(numWindows);
    
    for clustIdx=1:numClust
      
      if clustIdx<numClust
        locIdxDum=locIdx(1:tempperclust);
        locIdx=locIdx(tempperclust+1:end);
      else
        locIdxDum=locIdx;
      end
      locIdxDum=locIdxDum(~isnan(loc{T}(locIdxDum))); % remove non-windows (masked data/ NaNs)
      clustID(locIdxDum,T)=clustIdx;
      [trldum tIdxDum]=ind2sub([numTrl, winPerTrial],locIdxDum);
      clust{clustIdx,T}.linIdx=locIdxDum;
      clust{clustIdx,T}.trl=trldum;
      clust{clustIdx,T}.tIdx=tIdxDum;
      clust{clustIdx,T}.numWindows=numel(locIdxDum);
      
    end
  end
  
else %only reconstruct clust ID and initial cost
  clustID=nan(numWindows,nPT);
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

% initialize the cost function
if verbose
  fprintf('\nInitializing cost matrices...')
end
D=zeros(1,nPT);
N=winLen*prod(sz(3:end));
for T=1:nPT
  D(T)=0;
  for n=1:numClust
    N_c=clust{n,T}.numWindows;
    z_isum=nansum(z_i{T}(clust{n,T}.linIdx,:),1);
    if zscoreFlag
      clust{n,T}.varianceFac=N_c;
    else
      clust{n,T}.varianceFac=N_c*nanmean(nanmean(z_i{T}(clust{n,T}.linIdx,:).^2,2));
    end
    varianceFac=clust{n,T}.varianceFac;
    clust{n,T}.z_isum=z_isum;
    clust{n,T}.noNanCount=sum(~isnan(z_i{T}(clust{n,T}.linIdx,:)));
    if nanFlag
      % correct for NaNs
      nanFac=N_c./clust{n,T}.noNanCount;
      clust{n,T}.tempCost=(varianceFac*N_c/(N_c-1))-((nanFac.^2.*z_isum)*z_isum.')/((N-1)*(N_c-1));
    else
      clust{n,T}.tempCost=(varianceFac*N_c/(N_c-1))-(z_isum*z_isum.')/((N-1)*(N_c-1));
    end
    
    D(T)=D(T)+clust{n,T}.tempCost;
  end
end

if verbose
  fprintf('\nDone!\n')
end

% sort the locs such that lowest temperatures have lowest cost
[D,srtIdx]=sort(D);
loc=loc(srtIdx);
clust=clust(:,srtIdx);
clustID=clustID(:,srtIdx);
z_i=z_i(srtIdx);

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

swapCount=0;
tunnelCount=0;
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
  swapCount=swapCount+1;
  saveCount=saveCount+outputFlag;
  tunnelCount=tunnelCount+1;
  Tchangecount=Tchangecount+1;
  rejcount(2,:)=rejcount(2,:)+1;
  
  %determine to do cluster switch or window shift of sample
  pshift=.5;
  if numClust==1 %only one cluster means: never change clusters
    pshift=2;
  end
  
  for T=1:nPT
    linIdx=randval(numWindows);
    
    while maskFlag && isnan(loc{T}(linIdx)) % do not try to move pruned locs
      linIdx=randval(numWindows);
    end
    clustIdx=clustID(linIdx,T);
    
    shift=rand<pshift;
    
    if shift %shift a window (no cluster swap)
      [trl tIdx]=ind2sub([numTrl,winPerTrial],linIdx);
      
      pLoc=loc{T}(trl,tIdx);
      locChange=0;
      loopcount=0;
      while ~locChange && loopcount<11 %find a possible step; Brute force for maximally 10 times
        dir=((rand>0.5)-0.5)*2*round(stepSz*rand);
        nLoc=pLoc+dir;
        % preferably have it within data border, if not, at least move towards border
        leftLim= nLoc > 0 || nLoc > pLoc;
        rightLim= size(dat,2)-nLoc>=winLen  || nLoc<pLoc;
        locChange= leftLim && rightLim;
        
        if locChange && winPerTrial>1
          %check guard:
          otherWinSel=true(winPerTrial,1);
          otherWinSel(tIdx)=false;
          minDist=min(abs(loc{T}(trl,otherWinSel)-nLoc));
          if ~isnan(minDist) % if minDist is NaN this means there are no other windows in this trial
            locChange= minDist >= guard;
          end
        end
        
        % also check for distance from mask
        if locChange && maskFlag
          selDum=nLoc:min(nLoc+min(winLen,guard)-1,size(dat,2));
          selDum=selDum(selDum>0);
          maskDum=sum(mask(trl,selDum))<1;
          locChange=maskDum;
        end
        loopcount=loopcount+1;
      end
      
      
      
      if locChange %valid step is found, now evaluate change in cost function
        
        N_c=clust{clustIdx,T}.numWindows;
        %       pcost=D(T);
        pcost=clust{clustIdx,T}.tempCost;
        s_dum=nan(winLen,sz(3:end));
        if nLoc <1
          selVec=1:nLoc+winLen-1;
          selVecS=2-nLoc:winLen;
        elseif (nLoc + winLen -1)>sz(2)
          selVec=nLoc:sz(2);
          selVecS=1:numel(selVec);
        else % no problems
          selVec=nLoc:nLoc+winLen-1;
          selVecS=1:winLen;
        end
        s_dum(selVecS,:)=dat(trl,selVec,:);
        
        if zscoreFlag
          nZ=bsxfun(@times,kernel,(s_dum-nanmean(s_dum(:)))/nanstd(s_dum(:),1));
          varianceFac=N_c;
        else
          nZ=bsxfun(@times,kernel,s_dum);
          varianceFac=clust{clustIdx,T}.varianceFac;
          pZ=z_i{T}(linIdx,:);
          varianceFac=varianceFac-nanmean(pZ.^2)+nanmean(z_dum(:).^2);
        end
        
        
        pZ=z_i{T}(linIdx,:);
        if nanFlag
          noNanCountdum=clust{clustIdx,T}.noNanCount-isnan(pZ)+isnan(nZ(:)');
          pZ(isnan(pZ))=0;
          nZ(isnan(nZ))=0;
          z_sumdum=clust{clustIdx,T}.z_isum-pZ+nZ(:)';
          
          nanFac=N_c./noNanCountdum;
          ncost=(varianceFac*N_c/(N_c-1))-((nanFac.^2.*z_sumdum)*z_sumdum.')/((N-1)*(N_c-1));
        else
          z_sumdum=clust{clustIdx,T}.z_isum-pZ+nZ(:)';
          ncost=(varianceFac*N_c/(N_c-1))-(z_sumdum*z_sumdum.')/((N-1)*(N_c-1));
        end
        
        
        
        % change in cost function
        cVal=(pcost-ncost);
      else %could not find a valid window shift -> "reject", i.e. keep everything as-is
        cVal=-inf;
      end
      
      % Determine whether to accept or reject window shift
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
        loc{T}(trl,tIdx)=nLoc;
        clust{clustIdx,T}.z_isum=z_sumdum;
        clust{clustIdx,T}.varianceFac=varianceFac;
        if nanFlag
          clust{clustIdx,T}.noNanCount=noNanCountdum;
        end
        clust{clustIdx,T}.tempCost=ncost;
        D(T)=D(T)-cVal;
        z_i{T}(linIdx,:)=nZ(:);
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
      relinIdx=clust{clustIdx,T}.linIdx==linIdx;
      nclustIdx=randval(numClust);
      while nclustIdx == clustIdx
        nclustIdx=randval(numClust);
      end
      pcost=clust{clustIdx,T}.tempCost+clust{nclustIdx,T}.tempCost;
      
      N_c=[clust{clustIdx,T}.numWindows clust{nclustIdx,T}.numWindows]+[-1 1];
      
      cZ=z_i{T}(linIdx,:);
      if nanFlag
        
        cZ(isnan(cZ))=0;
        z_sumdum=[clust{clustIdx,T}.z_isum-cZ; clust{nclustIdx,T}.z_isum+cZ];
        
        noNanCountdum=[clust{clustIdx,T}.noNanCount-isnan(z_i{T}(linIdx,:)); clust{nclustIdx,T}.noNanCount+isnan(z_i{T}(linIdx,:))];
        nanFac=bsxfun(@rdivide,N_c.',noNanCountdum);
        z2dum=[(nanFac(1,:).^2.*z_sumdum(1,:))*z_sumdum(1,:).' (nanFac(2,:).^2.*z_sumdum(2,:))*z_sumdum(2,:).'];
      else
        z_sumdum=[clust{clustIdx,T}.z_isum-cZ; clust{nclustIdx,T}.z_isum+cZ];
        z2dum=[z_sumdum(1,:)*z_sumdum(1,:).' z_sumdum(2,:)*z_sumdum(2,:).'];
      end
      
      if zscoreFlag
        varianceFac=N_c;
      else
        varianceFac=[clust{clustIdx,T}.varianceFac clust{nclustIdx,T}.varianceFac];
        varianceFac=varianceFac+[-nanmean(cZ.^2) nanmean(cZ.^2)];
      end
      
      ncost=((varianceFac.*N_c./(N_c-1))-z2dum./((N-1)*(N_c-1)));
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
        clustID(linIdx,T)=nclustIdx;
        clust{clustIdx,T}.linIdx=clust{clustIdx,T}.linIdx(~relinIdx);
        clust{nclustIdx,T}.linIdx=[clust{nclustIdx,T}.linIdx linIdx];
        clust{clustIdx,T}.z_isum=z_sumdum(1,:);
        clust{nclustIdx,T}.z_isum=z_sumdum(2,:);
        clust{clustIdx,T}.tempCost=ncost(1);
        clust{nclustIdx,T}.tempCost=ncost(2);
        clust{clustIdx,T}.numWindows=clust{clustIdx,T}.numWindows-1;
        clust{nclustIdx,T}.numWindows=clust{nclustIdx,T}.numWindows+1;
        clust{nclustIdx,T}.trl=[clust{nclustIdx,T}.trl clust{clustIdx,T}.trl(relinIdx)];
        clust{nclustIdx,T}.tIdx=[clust{nclustIdx,T}.tIdx clust{clustIdx,T}.tIdx(relinIdx)];
        clust{clustIdx,T}.trl=clust{clustIdx,T}.trl(~relinIdx);
        clust{clustIdx,T}.tIdx=clust{clustIdx,T}.tIdx(~relinIdx);
        clust{clustIdx,T}.varianceFac=varianceFac(1);
        clust{nclustIdx,T}.varianceFac=varianceFac(2);
        D(T)=D(T)-cVal;
        
        if costMin(T)>D(T)
          if mincostTot>D(T);
            mincostTot=D(T);
            tloc=loc{T};
            tclustID=clustID(:,T);
            tclust=clust(:,T);
            for nn=1:numClust
              clustnumel(nn)=tclust{nn}.numWindows;
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
  if nPT>1 && swapCount >= TswapInterval;
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
    swapCount=0;
  end
  
  % Every couple of interations, try moving the entire FoV (moving all
  % windows at once)
  if tunnelCount >= FoVShiftInterval
    tunnelCount=0;
    for T=1:nPT
      % try 'tunneling' all clusters in a random order
      for clustIdx=randperm(numClust)
        linIdx=clust{clustIdx,T}.linIdx;
        stepFoV=sign(rand-.5)*randperm(ceil(guard/2),1);
        pLoc=loc{T};
        nLoc=pLoc;
        nLoc(linIdx)=nLoc(linIdx)+stepFoV;
        
        % check whether there are clashes with other clusters
        if size(nLoc,2)>1 && numClust>1
          winDist=diff(sort(nLoc,2),1,2);
          tolerance = 0; % no collissions are tolerable
          invalid=nanmean(winDist(:)<guard)>tolerance;
          if invalid
            continue % do not shift FoV
          end
        end
        
        locMask=true(size(pLoc));
        locMask(linIdx)=false;
        extractLoc=nLoc;
        extractLoc(locMask)=nan;
        
        % do NOT move windows outside of the data
        sel=extractLoc<1 | extractLoc+winLen-1>sz(2);
        extractLoc(sel)=pLoc(sel);
        
        % calculate new z_isum
        extractCfg=[];
        extractCfg.loc=extractLoc;
        extractCfg.winLen=winLen;
        extractCfg.numWindows=numWindows;
        
        % check whether shift is not clashing with the mask
        if maskFlag
          [maskdum]=bg_swm_extract(extractCfg,mask);
          maskPercent= nanmean(maskdum(:,:));
          tolerance=.1; % up to 10% mask is allowed
          if  any(maskPercent>tolerance)
            continue % do not shift FoV
          end
        end
        
        [s_N,z_N]=bg_swm_extract(extractCfg,dat);
        
        % cut out masked trials
        s_New=nan([numel(linIdx),winLen,sz(3:end)]);
        z_New=s_New;
        s_New(:,:)=s_N(linIdx,:);
        z_New(:,:)=z_N(linIdx,:);
        oldnanFlag=nanFlag;
        nanFlag=any(isnan(z_New(:)));
%         
%         if nanFlag
%           %check whether the shift of FoV is not too much. I.e. it should not
%           %move it from the data into nothingness
%           nanPercent= sort(mean(isnan(s_New(:,:))),'descend');
%           
%           %do not tolerate when 5% of windows contain at least 25% NaNs
%           %(at the same positions, indicating SWM is walking away from the data)
%           toleranceLen=ceil(.05*numel(nanPercent));
%           invalid=mean(nanPercent(1:toleranceLen))>.25;
%           if invalid
%             nanFlag=oldnanFlag;
%             continue % do not shift FoV
%           end
%         end
        
        
        
        N_c=clust{clustIdx,T}.numWindows;
        pcost=clust{clustIdx,T}.tempCost;
        
        if zscoreFlag
          z_dum=bsxfun(@times,kernel,z_New);
          varianceFac=N_c;
        else
          z_dum=bsxfun(@times,kernel,s_New);
          varianceFac=N_c*nanmean(nanmean(z_dum.^2,2));
        end
        
        nZ=z_dum;
        z_sumdum=sum(nZ(:,:));
        if nanFlag          
          nZ(isnan(nZ))=0;
          noNanCountdum=sum(~isnan(z_dum(:,:)),1);
          nanFac=N_c./noNanCountdum;
          z_sumdum=sum(nZ(:,:));
          z_sumdum=z_sumdum.*nanFac;
        end
        ncost=(varianceFac*N_c/(N_c-1))-(z_sumdum*z_sumdum.')/((N-1)*(N_c-1));
    
        
                
        % change in cost function
        cVal=(pcost-ncost);
        
        
        % Determine whether to accept or reject FoV shift
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
          loc{T}=nLoc;
          clust{clustIdx,T}.z_isum=z_sumdum;
          clust{clustIdx,T}.varianceFac=varianceFac;
          if nanFlag
            clust{clustIdx,T}.noNanCount=noNanCountdum;
          end
          clust{clustIdx,T}.tempCost=ncost;
          D(T)=D(T)-cVal;
          z_i{T}(linIdx,:)=z_dum(:,:);
          if costMin(T)>D(T)
            if mincostTot>D(T);
              mincostTot=D(T);
              tloc=loc{T};
            end
            costMin(T)=D(T);
            
          end
        else
          nanFlag=oldnanFlag;
        end
      end
    end
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
      if compatibilityFlag
        set(get(hleg,'title'),'string','Tfac')
      end
      plotLegend=0;
    end
    
    if multiDim
      for n=1:min(numClust,3)
        subplot(min(numClust,3),2,sub2ind([2 min(numClust,3)],2,n))
        imDum=datMean+datStd*bsxfun(@rdivide,reshape(clust{n,1}.z_isum/clust{n,1}.numWindows,winLen,[]),kernel)';
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
        subplot(1,2,2)
        for TT=1:numel(plotselT)
          plot(clust{n,datMean+datStd*plotselT(TT)}.z_isum/clust{n,plotselT(TT)}.numWindows./kernel,'color',colorOrderShape(TT,:),'linewidth',2)
          hold on
        end
        hold off
        hleg2=legend(num2str(Tfac(plotselT)','%1.2e'));
        if compatibilityFlag
          set(get(hleg2,'title'),'string','Tfac')
        end
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
            plot(datMean+datStd*clust{n,plotselT(TT)}.z_isum/clust{n,plotselT(TT)}.numWindows./kernel,'linewidth',2)
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
      
      locMask=true(size(tloc));
      locMask(tclust{n}.linIdx)=false;
      
      extractLoc=tloc;
      extractLoc(locMask)=nan;
      % calculate new z_isum
      extractCfg=[];
      extractCfg.best_loc=extractLoc;
      extractCfg.winLen=winLen;
      extractCfg.numWindows=numWindows;
      [s_Out,z_Out]=bg_swm_extract(extractCfg,dat);
      
      
      cfg.best_s(n,:)=reshape(nanmean(s_Out,1),[],1);
      cfg.best_z(n,:)=reshape(nanmean(z_Out,1),[],1);
      
      cfg.costDistr{n}=cost_i(z_score(s_Out));
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
      cleanFields={'cc', 'costFinal', 'costDistr', 'costCoM', 'cm','costCoM_i'};
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

function [loc, winPerTrial, stepSz]=initloc(guard,winLen,data,winPerTrial,mask)
% initialize window locations
len=size(data);
maxWin= floor((len(2)-winLen+1)/(guard+1)*.75); % 75% data coverage; (100% would be too dense, no movement for the windows)
if nargin<4 || isempty(winPerTrial)
  winPerTrial=maxWin;
else
  if winPerTrial>maxWin
    winPerTrial=maxWin;
  end
end

if nargin<5 || isempty(mask)
  maskFlag=false;
else
  maskFlag=true;
end

if maskFlag % cut off leading and trailing mask
  startIdx=find(mean(mask)<.1,1,'first');
  endIdx=find(mean(mean(mask(:,:,:)),3)<.1,1,'last');
  len(2)=endIdx-startIdx+1;
end
  

loc=nan(size(data,1),winPerTrial);

nEmptySpace=len(2)-winLen-(winPerTrial-1)*guard;
guards=[0:winPerTrial-1]*guard;

for n=1:len(1)
  empty=sort(randperm(nEmptySpace,winPerTrial));
  winStarts=empty+guards;  
  loc(n,:)=sort(winStarts);  
end

if maskFlag
  loc=loc+startIdx-1;
end

stepSz=ceil(nEmptySpace/winPerTrial);

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
