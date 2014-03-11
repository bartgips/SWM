function [z, s, cfg]=bg_swm_extract(cfg,dat)
% [z, s, cfg]=bg_swm_extract(dat,cfg)
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
% last edit: 6 March 2014

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


if isfield(cfg,'newLen')
  newLen=cfg.newLen;
else
  newLen=cfg.winLen;
end

if isfield(cfg,'best_loc')
  loc=cfg.best_loc;
else loc=cfg.loc;
end
winLen=cfg.winLen;

if isfield(cfg,'best_clust')
  cfg.clust=cfg.best_clust;
  cfg.numclust=numel(cfg.clust);
end

if isfield(cfg,'clust')
  clustering=true;
  s=cell(1,cfg.numclust);
  for n=1:numel(s)
    s{n}=nan(cfg.clust{n}.numTemplates,newLen);
  end
else
  clustering=false;
  s=nan(cfg.numTemplates,newLen);
end
z=s;
addlen=round(newLen-winLen)/2;
cfg.addlen=addlen;
newLen=winLen+2*addlen;


if clustering
  for n=1:numel(s)
    for k=1:cfg.clust{n}.numTemplates
      trl=cfg.clust{n}.trl(k);
      tidx=cfg.clust{n}.tidx(k);
      if loc(trl,tidx)+winLen+addlen<=size(dat,2) && loc(trl,tidx)-addlen>0
        s{n}(k,:)=dat(trl,loc(trl,tidx)-addlen:loc(trl,tidx)+winLen+addlen-1);
      elseif loc(trl,tidx)-addlen>0
        num=size(dat,2)-loc(trl,tidx)+1;
        s{n}(k,1:num)=dat(trl,loc(trl,tidx):end);
      elseif loc(trl,tidx)+addlen<=size(dat,2)
        num=1-(loc(trl,tidx)-addlen);
        s{n}(k,num+1:end)=dat(trl,1:loc(trl,tidx)+winLen+addlen-1);
      end
    end
    z{n}=zscore(s{n},1,2);
  end
else
  for k=1:cfg.numTemplates
    [trl, tidx]=ind2sub(size(loc(:,1:end-1)),k);
    if loc(trl,tidx)+winLen+addlen<=size(dat,2) && loc(trl,tidx)-addlen>0
      s(k,:)=dat(trl,loc(trl,tidx)-addlen:loc(trl,tidx)+winLen+addlen-1);
    elseif loc(trl,tidx)-addlen>0
      num=size(dat,2)-loc(trl,tidx)+1;
      s(k,1:num)=dat(trl,loc(trl,tidx):end);
    elseif loc(trl,tidx)+addlen<=size(dat,2)
      num=1-(loc(trl,tidx)-addlen);
      s(k,num+1:end)=dat(trl,1:loc(trl,tidx)+winLen+addlen-1);
    end
  end
  z=zscore(s,1,2);
end


function z=zscore(x,flag,dim)
% 
if nargin<3
  dim=1;
end
if nargin<2
  flag=0;
end

sdx=std(x,flag,dim);
mux=mean(x,dim);

z=bsxfun(@minus,x,mux);
z=bsxfun(@rdivide,z,sdx);
