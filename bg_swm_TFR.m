function [TFR, f, sd, sem, PLV]=bg_swm_TFR(cfg, freqoi, dat)
% [TFR, f, sd, sem, PLV]=bg_extract_TFR(cfg, freqoi, dat)
%
% input:
% cfg:    output from bg_SWM
% freqoi: frequencies of interest
% dat:    optional; data that the cfg belongs to. (only relevant if cfg
%         does not contain .varName and .fName
% 
% last edit: 6 March 2014


if nargin<3
  dum=load(cfg.fName, cfg.varName);
  eval(['dat=dum.' cfg.varName ';'])
end

if iscolumn(dat)
  dat=dat.';
end

winLen=cfg.winLen;

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

if isfield(cfg,'best_loc')
  loc=cfg.best_loc;
else
  loc=cfg.loc;
end

if isfield(cfg,'newLen')
  newLen=cfg.newLen;
  addLen=round((newLen-winLen)/2);
  cfg.addLen=addLen;
  newLen=winLen+2*addLen;
  cfg.newLen=newLen;
  winLen=newLen;
else
  addLen=0;
end


if isfield(cfg,'best_clust')
  cfg.clust=cfg.best_clust;
elseif ~isfield(cfg,'clust')
  cfg.numClust=1;
  cfg.clust{1}.linIdx=1:numel(loc);
  [trl, tIdx]=ind2sub(size(loc),1:numel(loc));
  cfg.clust{1}.trl=trl;
  cfg.clust{1}.tIdx=tIdx;
end

if isfield(cfg,'derivative')
  derivative=cfg.derivative;
else
  derivative=0;
end

if derivative
  dat=diff(dat,1,2);
end



t=[1:size(dat,2)]/cfg.fs;


% freqoi=[1:2:100];
timoi=t;
groupsz=min(50,size(dat,1));
tfr=nan(groupsz,numel(freqoi),numel(timoi));
numgroups=ceil(size(dat,1)/groupsz);
if verbose
  fprintf('\nGenerating TFR of data in groups of %d trials/channels...\n',groupsz)
  reverseStr='';
end
TFR=zeros(numel(freqoi),winLen,cfg.numClust);
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
  for kk=1:cfg.numClust
    sel=find((cfg.clust{kk}.trl<(n*groupsz+1) & cfg.clust{kk}.trl>(n-1)*groupsz));
    trl=cfg.clust{kk}.trl(sel)-(n-1)*groupsz;
    tIdx=cfg.clust{kk}.tIdx(sel);
    
    % remove NaN locs
    lidxdum=sub2ind(size(loc),trl+(n-1)*groupsz,tIdx);
    nansel=isnan(loc(lidxdum));
    sel=sel(~nansel);
    trl=cfg.clust{kk}.trl(sel)-(n-1)*groupsz;
    tIdx=cfg.clust{kk}.tIdx(sel);
    for k=1:numel(sel);
      strt=loc(trl(k)+(n-1)*groupsz,tIdx(k))-addLen;
      fin=strt+winLen-1;
      if strt<1
        num=1-strt;
        avgdum(:,num+1:end,kk)=avgdum(:,num+1:end,kk)+row(squeeze(tfrnan(trl(k),:,1:fin)));
        sigdumC=0*TFRdum;
        sigdumC(:,num+1:end,kk)=row(squeeze(tfr(trl(k),:,1:fin)));
        sigdum=abs(sigdumC).^2;
        if nargout>2 %calculate running std: Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd ed; http://www.johndcook.com/standard_deviation.html
          if n==1 && k==1 %initialize M1=x1;
            stddumM(:,:,kk)=sigdum(:,:,kk);
          end
          stddumN(:,num+1:end,kk)=stddumN(:,num+1:end,kk)+row(squeeze(tfrnan(trl(k),:,1:fin)));
          stddumMp=stddumM(:,:,kk);
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumM(stddumN==0)=0; % correct for cases where no data is added yet (otherwise NaNs show up)
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumMp(:,:,kk)).*(sigdum(:,:,kk)-stddumM);
          
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
          if n==1 && k==1 %initialize M1=x1;
            stddumM(:,:,kk)=sigdum(:,:,kk);
          end
          stddumN(:,1:num,kk)=stddumN(:,1:num,kk)+row(squeeze(tfrnan(trl(k),:,strt:size(tfr,3))));
          stddumMp=stddumM(:,:,kk);
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumM(stddumN==0)=0; % correct for cases where no data is added yet (otherwise NaNs show up)
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumMp(:,:,kk)).*(sigdum(:,:,kk)-stddumM);
          
        end
        
      else
        sigdumC=0*TFRdum;
        sigdumC(:,:,kk)=row(squeeze(tfr(trl(k),:,strt:fin)));
        sigdum=abs(sigdumC).^2;
        
        
        TFRdum(:,:,kk)=TFRdum(:,:,kk)+sigdum(:,:,kk);
        TFRCdum(:,:,kk)=TFRCdum(:,:,kk)+sigdumC(:,:,kk);
        TFRAdum(:,:,kk)=TFRAdum(:,:,kk)+abs(sigdumC(:,:,kk));
        avgdum(:,:,kk)=avgdum(:,:,kk)+row(squeeze(tfrnan(trl(k),:,strt:fin)));
        
        if nargout>2 %calculate running std
          if n==1 && k==1 %initialize M1=x1;
            stddumM(:,:,kk)=sigdum(:,:,kk);
          end
          stddumN(:,:,kk)=stddumN(:,:,kk)+row(squeeze(tfrnan(trl(k),:,strt:fin)));
          stddumMp=stddumM(:,:,kk);
          stddumM(:,:,kk)=stddumM(:,:,kk)+(sigdum(:,:,kk)-stddumM(:,:,kk))./stddumN(:,:,kk);
          stddumM(stddumN==0)=0; % correct for cases where no data is added yet (otherwise NaNs show up)
          stddumS(:,:,kk)=stddumS(:,:,kk)+(sigdum(:,:,kk)-stddumMp(:,:,kk)).*(sigdum(:,:,kk)-stddumM);
          
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
