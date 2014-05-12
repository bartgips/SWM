function [z, s, cfg]=bg_swm_extract(cfg, dat)
% [z, s, cfg]=bg_swm_extract(cfg, dat)
%
% input:
% cfg:    output from bg_SWM
%         can contain optional new field: .newLen;
%         this is the length of the windows that should be cut out of the
%         dataset (symetrically around the orignal window positions)
%
% dat:    optional; data that the cfg belongs to. (only relevant if cfg
%         does not contain .varname and .fname
%
%
%
% last edit: March 2014

if nargin<2
  dum=load(cfg.fname, cfg.varname);
  eval(['dat=dum.' cfg.varname ';'])
else
  cfg.fname='function input';
  cfg.varname='N/A';
end

if iscolumn(dat)
  dat=dat.';
end


%% preprocessing (filtering)
nanSel=isnan(dat);
nanFlag=any(nanSel(:));

if any(isfield(cfg,{'Fbp','Fbs','Fhp','Flp'}))
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

if isfield(cfg,'best_loc')
  loc=cfg.best_loc;
else
  loc=cfg.loc;
end

if isfield(cfg,'best_clust')
  cfg.clust=cfg.best_clust;
  cfg.numclust=numel(cfg.clust);
end


if isfield(cfg,'clust') && numel(cfg.clust)>1
  clustering=true;
  s=cell(1,cfg.numclust);
  for n=1:numel(s)
    s{n}=nan(cfg.clust{n}.numTemplates,newLen);
  end
else
  clustering=false;
  s=nan(cfg.numTemplates,newLen);
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
  s=nan(numel(loc),newLen);
  
else
  manualLoc=false;
end


z=s;




if manualLoc
  for k=1:numel(loc);
    trl=real(loc(k));
    tidx=imag(loc(k));
    if tidx+winLen+addLen<=size(dat,2) && tidx-addLen>0
      s(k,:)=dat(trl,tidx-addLen:tidx+winLen+addLen-1);
    elseif tidx-addLen>0
      num=size(dat,2)-tidx+1;
      s(k,1:num)=dat(trl,tidx:end);
    elseif tidx+addLen+winLen<=size(dat,2)
      num=1-(tidx-addLen);
      s(k,num+1:end)=dat(trl,1:tidx+winLen+addLen-1);
      %     else
      %       strt=tidx+addLen;
      %       finish=size(dat,2)+tidx+addLen;
      %       s(k,strt+1:finish)=dat(trl,1:end);
    end
  end
  z=zscore(s,1,2);
  
else
  if clustering
    for n=1:numel(s)
      for k=1:cfg.clust{n}.numTemplates
        trl=cfg.clust{n}.trl(k);
        tidx=cfg.clust{n}.tidx(k);
        if loc(trl,tidx)+winLen+addLen<=size(dat,2) && loc(trl,tidx)-addLen>0
          s{n}(k,:)=dat(trl,loc(trl,tidx)-addLen:loc(trl,tidx)+winLen+addLen-1);
        elseif loc(trl,tidx)-addLen>0
          num=size(dat,2)-loc(trl,tidx)+1;
          s{n}(k,1:num)=dat(trl,loc(trl,tidx):end);
        elseif loc(trl,tidx)+addLen<=size(dat,2)
          num=1-(loc(trl,tidx)-addLen);
          s{n}(k,num+1:end)=dat(trl,1:loc(trl,tidx)+winLen+addLen-1);
        end
      end
      z{n}=zscore(s{n},1,2);
    end
  else
    for k=1:cfg.numTemplates
      [trl, tidx]=ind2sub(size(loc(:,1:end-1)),k);
      if loc(trl,tidx)+winLen+addLen<=size(dat,2) && loc(trl,tidx)-addLen>0
        s(k,:)=dat(trl,loc(trl,tidx)-addLen:loc(trl,tidx)+winLen+addLen-1);
      elseif loc(trl,tidx)-addLen>0
        num=size(dat,2)-loc(trl,tidx)+1;
        s(k,1:num)=dat(trl,loc(trl,tidx):end);
      elseif loc(trl,tidx)+addLen+winLen<=size(dat,2)
        num=1-(loc(trl,tidx)-addLen);
        s(k,num+1:end)=dat(trl,1:loc(trl,tidx)+winLen+addLen-1);
        %     else
        %       strt=loc(trl,tidx)+addLen;
        %       finish=size(dat,2)+loc(trl,tidx)+addLen;
        %       s(k,strt+1:finish)=dat(trl,1:end);
      end
    end
    z=zscore(s,1,2);
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
