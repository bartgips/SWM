function [s, z, cfg]=bg_swm_extract(cfg, dat)
% [s, z, cfg]=bg_swm_extract(cfg, dat)
%
% input:
% cfg:    output from bg_SWM
%         can contain optional new field: .newLen;
%         this is the length of the windows that should be cut out of the
%         dataset (symetrically around the orignal window positions)
%
% dat:    optional; data that the cfg belongs to. (only relevant if cfg
%         does not contain .varName and .fName
%
%
%
% last edit: March 2014

%% check validity of input fields
validInp={'best_s';'best_z';'best_clust';'best_clustID';'best_loc';'clust';...
  'costDistr';'costMin';'costFinal';'costTotal';...
  'costTotal_end';'costTotal_undSamp';'debug';'dispPlot';'Fbs';...
  'Fbp';'Fhp';'FhpFac';'Flp';'fName';'fs';'fullOutput';'guard';'kernel';...
  'konstant';'loc';'mask';'normalize';'nPT';'numClust';'numIt';'numIter';'numWindows';...
  'winPerTrial';'outputFile';'stepSz';'Tfac';'varName';'verbose';'winLen';'winLenFac';'zscore'};

% fix case of inputs
try
  [~,cfg]=isfieldi(cfg,validInp);
catch
  error('cfg contains conflicting fields')
end





if nargin<2
  dum=load(cfg.fName, cfg.varName);
  eval(['dat=dum.' cfg.varName ';'])
else
  cfg.fName='function input';
  cfg.varName='N/A';
end

if iscolumn(dat)
  dat=dat.';
end

sz=size(dat);
if numel(sz)<3
  sz(3)=1;
end


%% preprocessing (filtering)
nanSel=isnan(dat);
nanFlag=any(nanSel(:));

[val,cfgOut]=isfieldi(cfg,{'Fbp','Fbs','Fhp','Flp'});
if any(val)
  cfg=cfgOut;
  disp('Applying frequency filters...')
  %   filttype='fir';
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
  if cfg.Fhp>0
  Fhp=cfg.Fhp;
  Fhp=Fhp(:);
  for freq=1:size(Fhp,1)
    for n=1:prod(sz(3:end))
    dat(:,:,n)=ft_preproc_highpassfilter(dat(:,:,n), fs, Fhp(freq),[],filttype);
    end
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

dat(nanSel)=nan;
%%
if isfield(cfg,'newLen')
  newLen=cfg.newLen;
else
  newLen=cfg.winLen;
end

winLen=cfg.winLen;
addLen=round((newLen-winLen)/2);
cfg.addLen=addLen;
newLen=winLen+2*addLen;
cfg.newLen=newLen;

if isfield(cfg,'best_loc') && ~all(isnan(cfg.best_loc(:)))
  loc=cfg.best_loc;
else
  loc=cfg.loc;
  if iscell(loc)
    loc=loc{1};
  end
end

if isfield(cfg,'best_clust')
  cfg.clust=cfg.best_clust;
  cfg.numclust=numel(cfg.clust);
end

validClustFields={'linIdx','trl','tIdx'};
if isfield(cfg,'clust') && numel(cfg.clust)>1
  clustering=true;
  s=cell(1,cfg.numclust);
  for n=1:numel(s)
    s{n}=nan([cfg.clust{n}.numWindows,newLen,sz(3:end)]);
  end
  for n=1:numel(cfg.clust)
    [~,cfg.clust{n}]=isfieldi(cfg.clust{n},validClustFields);
  end
else
  clustering=false;
  s=nan([cfg.numWindows,newLen,sz(3:end)]);
end



if isfield(cfg,'loc_manual') % manual input of window locations.
  % loc is now a vector of complex values. Real part: trial index, Imag
  % part: time index.
  loc=cfg.loc_manual;
  if ~iscomplex(loc)
    error('loc should be a vector of complex values. Real part: trial index, Imag: time index')
  end
  % part: time index.
  manualLoc=true;
  clustering=false;
  s=nan([numel(loc),newLen,sz(3:end)]);
  
else
  manualLoc=false;
end

if ~isfield(cfg,'numWindows')
  cfg.numWindows=numel(loc);
end

if nargout>1
  z=s;
end




if manualLoc
  for k=1:numel(loc);
    trl=real(loc(k));
    tIdx=imag(loc(k));
    if tIdx+winLen+addLen<=size(dat,2) && tIdx-addLen>0
      s(k,:)=reshape(dat(trl,tIdx-addLen:tIdx+winLen+addLen-1,:),[],1);
    elseif tIdx-addLen>0
      num=size(dat,2)-tIdx+1;
      s(k,1:num,:)=dat(trl,tIdx:end,:);
    elseif tIdx+addLen+winLen<=size(dat,2)
      num=1-(tIdx-addLen);
      s(k,num+1:end,:)=dat(trl,1:tIdx+winLen+addLen-1,:);
      %     else
      %       strt=tIdx+addLen;
      %       finish=size(dat,2)+tIdx+addLen;
      %       s(k,strt+1:finish)=dat(trl,1:end);
    end
    if nargout>1
      z(k,:)=(s(k,:)-nanmean(s(k,:)))/nanstd(s(k,:));
    end
  end
  
else
  if clustering
    for n=1:numel(s)
      for k=1:cfg.clust{n}.numWindows
        trl=cfg.clust{n}.trl(k);
        tIdx=cfg.clust{n}.tIdx(k);
        if loc(trl,tIdx)+winLen+addLen<=size(dat,2) && loc(trl,tIdx)-addLen>0
          s{n}(k,:)=reshape(dat(trl,loc(trl,tIdx)-addLen:loc(trl,tIdx)+winLen+addLen-1,:),[],1);
        elseif loc(trl,tIdx)-addLen>0
          num=size(dat,2)-loc(trl,tIdx)+1;
          selVec=false(size(s{n}));
          selVec(k,1:num,:)=true;
          s{n}(selVec)=reshape(dat(trl,loc(trl,tIdx):end,:),[],1);
        elseif loc(trl,tIdx)+addLen<=size(dat,2)
          num=1-(loc(trl,tIdx)-addLen);
          selVec=false(size(s{n}));
          selVec(k,num+1:end,:)=true;
          s{n}(selVec)=reshape(dat(trl,1:loc(trl,tIdx)+winLen+addLen-1,:),[],1);
        end
      end
      if nargout>1
        z{n}=zscore(s{n},0,2);
      end
    end
  else
    for k=1:cfg.numWindows
      [trl, tIdx]=ind2sub(size(loc),k);
      if loc(trl,tIdx)+winLen+addLen<=size(dat,2) && loc(trl,tIdx)-addLen>0
        s(k,:)=reshape(dat(trl,loc(trl,tIdx)-addLen:loc(trl,tIdx)+winLen+addLen-1,:),[],1);
      elseif loc(trl,tIdx)-addLen>0 % start of window fits
        num=size(dat,2)-loc(trl,tIdx)+1+addLen;
        s(k,1:num,:)=dat(trl,loc(trl,tIdx)-addLen:end,:);
      elseif loc(trl,tIdx)+addLen+winLen<=size(dat,2) %end of window fits
        num=1-(loc(trl,tIdx)-addLen);
        s(k,num+1:end,:)=dat(trl,1:loc(trl,tIdx)+winLen+addLen-1,:);
        %     else
        %       strt=loc(trl,tIdx)+addLen;
        %       finish=size(dat,2)+loc(trl,tIdx)+addLen;
        %       s(k,strt+1:finish)=dat(trl,1:end);
      end
      if nargout>1
        z(k,:)=(s(k,:)-nanmean(s(k,:)))/nanstd(s(k,:));
      end
    end
  end
end


function z=zscore(x,flag,dim)
%
if nargin<3
  dim=1;
end
if nargin<2
  flag=0;
end

sdx=nanstd(x,flag,dim);
mux=nanmean(x,dim);

z=bsxfun(@minus,x,mux);
z=bsxfun(@rdivide,z,sdx);