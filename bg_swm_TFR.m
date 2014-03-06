function [TFR, f, sd, sem, PLV]=bg_swm_TFR(cfg, freqoi, dat)
% [TFR, f, sd, sem, PLV]=bg_extract_TFR(cfg, freqoi, dat)
%
% input:
% cfg:    output from bg_SWM
% freqoi: frequencies of interest
% dat:    optional; data that the cfg belongs to. (only relevant if cfg
%         does not contain .varname and .fname
% 
% last edit: 6 March 2014


if nargin<3
  dum=load(cfg.fname, cfg.varname);
  eval(['dat=dum.' cfg.varname ';'])
end

if iscolumn(dat)
  dat=dat.';
end

if isfield(cfg,'width')
  width=cfg.width;
else
  width=7;
end

if isfield(cfg,'verbose')
  verbose=cfg.verbose;
else
  verbose=1;
end

if isfield(cfg,'bestloc')
  loc=cfg.bestloc;
else
  loc=cfg.loc;
end

if isfield(cfg,'newlen') && isfield(cfg,'addlen')
  loc=loc-cfg.addlen;
  cfg.fitlen=cfg.newlen;
end

if isfield(cfg,'bestclust')
  cfg.clust=cfg.bestclust;
end

fitlen=cfg.fitlen;

t=[1:size(dat,2)]/cfg.fs;
tim=t(2:2:end);

% freqoi=[1:2:100];
timoi=t;
groupsz=min(100,size(dat,1));
tfr=nan(groupsz,numel(freqoi),numel(timoi));
numgroups=ceil(size(dat,1)/groupsz);
if verbose
  fprintf('\nGenerating TFR of data in groups of %d trials/channels...\n',groupsz)
  reverseStr='';
end
TFR=zeros(numel(freqoi),cfg.fitlen,cfg.numclust);
avg=TFR;
avgdum=TFR;
stddumM=TFR;
stddumS=stddumM;
stddumN=stddumM;
TFRCdum=0*TFR;
TFRAdum=0*TFR;

for n=1:numgroups
  if verbose
    msg=sprintf(' %d/%d\n', [n numgroups]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  if n==numgroups
    [dum, f, tim]=ft_specest_tfr(dat(((n-1)*groupsz+1):end,:), t, 'timeoi', timoi, 'freqoi', freqoi,'verbose',0,'width',width);
    tfr=nan*tfr;
    tfr(1:size(dum,1),1:numel(f),:)=dum;
  else
    [dum, f, tim]=ft_specest_tfr(dat([1:groupsz]+(n-1)*groupsz,:), t, 'timeoi', timoi, 'freqoi', freqoi, 'verbose', 0, 'width',width);
    
    tfr(1:groupsz,1:numel(f),:)=dum;
  end
  if n==1;
    TFR=TFR(1:numel(f),:,:);
  end
  tfrnan=~isnan(tfr);
  tfr(isnan(tfr))=0;
  TFRdum=0*TFR;
 
  
  avgdum=0*avgdum;
  for kk=1:cfg.numclust
    sel=find((cfg.clust{kk}.trl<(n*groupsz+1) & cfg.clust{kk}.trl>(n-1)*groupsz));
    trl=cfg.clust{kk}.trl(sel)-(n-1)*groupsz;
    tidx=cfg.clust{kk}.tidx(sel);
    for k=1:numel(sel);
      strt=loc(trl(k)+(n-1)*groupsz,tidx(k));
      fin=loc(trl(k)+(n-1)*groupsz,tidx(k))+fitlen-1;
      if strt<1
        num=1-strt;
        avgdum(:,num+1:end,kk)=avgdum(:,num+1:end,kk)+row(squeeze(tfrnan(trl(k),:,1:fin)));
        sigdumC=0*TFRdum;
        sigdumC(:,num+1:end,kk)=row(squeeze(tfr(trl(k),:,1:fin)));
        sigdum=abs(sigdumC).^2;
        if nargout>2 %calculate running std: Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd ed; http://www.johndcook.com/standard_deviation.html
          stddumN(:,num+1:end,kk)=stddumN(:,num+1:end,kk)+row(squeeze(tfrnan(trl(k),:,1:fin)));
          stddumMp=stddumM(:,:,kk);
          
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk)).*(sigdum(:,:,kk)-stddumMp);
          
        end
        
        
        TFRdum(:,num+1:end,kk)=TFRdum(:,num+1:end,kk)+sigdum(:,num+1:end,kk);
        TFRCdum(:,num+1:end,kk)=TFRCdum(:,num+1:end,kk)+sigdumC(:,num+1:end,kk);
        TFRAdum(:,num+1:end,kk)=TFRAdum(:,num+1:end,kk)+abs(sigdumC(:,num+1:end,kk));
      elseif fin> size(tfr,3)
        num=size(tfr,3)-strt+1;
        
        sigdumC=0*TFRdum;
        sigdumC(:,1:num,kk)=row(squeeze(tfr(trl(k),:,strt:size(tfr,3))));
        sigdum=abs(sigdumC).^2;
         
        
        TFRdum(:,1:num,kk)=TFRdum(:,1:num,kk)+sigdum(:,1:num,kk);
        TFRCdum(:,1:num,kk)=TFRCdum(:,1:num,kk)+sigdumC(:,1:num,kk);
        TFRAdum(:,1:num,kk)=TFRAdum(:,1:num,kk)+abs(sigdumC(:,1:num,kk));
        avgdum(:,1:num,kk)=avgdum(:,1:num,kk)+row(squeeze(tfrnan(trl(k),:,strt:size(tfr,3))));
        
        if nargout>2 %calculate running std
          stddumN(:,1:num,kk)=stddumN(:,1:num,kk)+row(squeeze(tfrnan(trl(k),:,strt:size(tfr,3))));
          stddumMp=stddumM(:,:,kk);
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk)).*(sigdum(:,:,kk)-stddumMp);
          
        end
        
      else
        sigdumC=0*TFRdum;
        sigdumC(:,:,kk)=row(squeeze(tfr(trl(k),:,loc(trl(k)+(n-1)*groupsz,tidx(k)):(loc(trl(k)+(n-1)*groupsz,tidx(k))+fitlen-1))));
        sigdum=abs(sigdumC).^2;
        
        
        TFRdum(:,:,kk)=TFRdum(:,:,kk)+sigdum(:,:,kk);
        TFRCdum(:,:,kk)=TFRCdum(:,:,kk)+sigdumC(:,:,kk);
        TFRAdum(:,:,kk)=TFRAdum(:,:,kk)+abs(sigdumC(:,:,kk));
        avgdum(:,:,kk)=avgdum(:,:,kk)+row(squeeze(tfrnan(trl(k),:,loc(trl(k)+(n-1)*groupsz,tidx(k)):(loc(trl(k)+(n-1)*groupsz,tidx(k))+fitlen-1))));
        
        if nargout>2 %calculate running std
          stddumN(:,:,kk)=stddumN(:,:,kk)+row(squeeze(tfrnan(trl(k),:,loc(trl(k)+(n-1)*groupsz,tidx(k)):(loc(trl(k)+(n-1)*groupsz,tidx(k))+fitlen-1))));
          stddumMp=stddumM(:,:,kk);
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk)).*(sigdum(:,:,kk)-stddumMp);
          
        end
        
      end
    end
    
  end
%   TFRdum=TFRdum./avgdum;
  TFRdum(isnan(TFRdum))=0; %nans occur when no values have been added in the average.
  TFR=TFR+TFRdum;
  avg=avg+avgdum;
  
end
PLV=(TFRCdum)./TFRAdum;
TFR=TFR./avg;

if nargout>2
  sd=sqrt(stddumS./(stddumN-1));
  if nargout>3
  sem=sd./sqrt(stddumN);
  end
end


function x=row(x)
if iscolumn(x)
  x=x.';
end
