function [z, s, cfg]=extract_shapes_symm(dat,cfg)
% [z, s, cfg]=extract_shapes_symm(dat,cfg)

if nargin<2
  cfg=dat;
  dum=load(cfg.fname, cfg.varname);
  eval(['dat=dum.' cfg.varname ';'])
else
  cfg.fname='function input';
  cfg.varname='N/A';
end

if iscolumn(dat)
  dat=dat.';
end


if isfield(cfg,'newlen')
  newlen=cfg.newlen;
else
  newlen=cfg.fitlen;
end

if isfield(cfg,'bestloc')
  loc=cfg.bestloc;
else loc=cfg.loc;
end
fitlen=cfg.fitlen;

if isfield(cfg,'bestclust')
  cfg.clust=cfg.bestclust;
  cfg.numclust=numel(cfg.clust);
end

if isfield(cfg,'clust')
  clustering=true;
  s=cell(1,cfg.numclust);
  for n=1:numel(s)
    s{n}=nan(cfg.clust{n}.numtemplates,newlen);
  end
else
  clustering=false;
  s=nan(cfg.numtemplates,newlen);
end
z=s;
addlen=round(newlen-fitlen)/2;
cfg.addlen=addlen;
newlen=fitlen+2*addlen;


if clustering
  for n=1:numel(s)
    for k=1:cfg.clust{n}.numtemplates
      trl=cfg.clust{n}.trl(k);
      tidx=cfg.clust{n}.tidx(k);
      if loc(trl,tidx)+fitlen+addlen<=size(dat,2) && loc(trl,tidx)-addlen>0
        s{n}(k,:)=dat(trl,loc(trl,tidx)-addlen:loc(trl,tidx)+fitlen+addlen-1);
      elseif loc(trl,tidx)-addlen>0
        num=size(dat,2)-loc(trl,tidx)+1;
        s{n}(k,1:num)=dat(trl,loc(trl,tidx):end);
      elseif loc(trl,tidx)+addlen<=size(dat,2)
        num=1-(loc(trl,tidx)-addlen);
        s{n}(k,num+1:end)=dat(trl,1:loc(trl,tidx)+fitlen+addlen-1);
      end
    end
    z{n}=zscore(s{n},1,2);
  end
else
  for k=1:cfg.numtemplates
    [trl, tidx]=ind2sub(size(loc(:,1:end-1)),k);
    if loc(trl,tidx)+fitlen+addlen<=size(dat,2) && loc(trl,tidx)-addlen>0
      s(k,:)=dat(trl,loc(trl,tidx)-addlen:loc(trl,tidx)+fitlen+addlen-1);
    elseif loc(trl,tidx)-addlen>0
      num=size(dat,2)-loc(trl,tidx)+1;
      s(k,1:num)=dat(trl,loc(trl,tidx):end);
    elseif loc(trl,tidx)+addlen<=size(dat,2)
      num=1-(loc(trl,tidx)-addlen);
      s(k,num+1:end)=dat(trl,1:loc(trl,tidx)+fitlen+addlen-1);
    end
  end
  z=zscore(s,1,2);
end