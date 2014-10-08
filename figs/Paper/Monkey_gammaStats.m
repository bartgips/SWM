clear all

addpath('/home/mrphys/bargip/OrientationDataMonkeys/matlab_preproc')


basedir={'/home/mrphys/bargip/OrientationDataMonkeys/Torque/Okkie';...
  '/home/mrphys/bargip/OrientationDataMonkeys/Torque/Spock'};

penSel={'009','142'};

Monkey={'Monkey O','Monkey S'};

chansel={1:16, 17:32};
% chansel={1:16, 1:16};

alignment=  [9 2 3 3; 141 0 2 2; 142 0 NaN 1; 143 1 7 NaN];
%% okkie Chansel
% 
% if strcmp(basedir(end-4:end),'Okkie')
%   chansel={1:16
%     1:48
%     1:48
%     1:48
%     1:16};
%   %1:16};
%   Fhp=25;
%   
%   Monkey='Okkie';
%   savdirbase='/home/mrphys/bargip/OrientationDataMonkeys/Torque/Okkie';
%   alignment=  [7 3 2 3
%     9 2 3 3
%     38 2 3 5
%     53 1 2 4];
% elseif strcmp(basedir(end-4:end),'Spock')
%   chansel={17:32
%     17:32
%     17:32
%     1:16
%     1:16
%     1:16
%     1:48
%     1:48
%     1:48
%     1:48};
%   Fhp=24;
%   Monkey='Spock';
%   
%   savdirbase='/home/mrphys/bargip/OrientationDataMonkeys/Torque/Spock';
%   alignment=  [160 5 4 4
%     161 5 4 3
%     162 3 3 3
%     139 0 NaN 8
%     141 0 2 2
%     142 0 NaN 1
%     143 1 7 NaN
%     144 2 5 5];
% end


%%
foi=[5:1:65];
[skew,r,Shape,t, mExtr, xCorrS, xCorrL]=deal(cell(numel(basedir),1));
aligOff=nan(2,3);
fnameBase='cfg_elec_';
shapeTim=cell(numel(basedir),1);
mainStats=cell(2,numel(basedir));
ttest=nan(1,numel(basedir));
[mu sem]=deal(nan(2,numel(basedir)));
%%

for monkeyIt=1:numel(basedir)
  pendir=dir(basedir{monkeyIt});
  pendir={pendir.name};
  pendir=bg_selectCell(pendir,'pen','~.mat');
  pendir=bg_selectCell(pendir,penSel{monkeyIt});
  savdir=fullfile(basedir{monkeyIt},pendir{1});
  
  if ~exist(savdir,'file')
    error(['Directory ' savdir ' does not exist'])
  end
  
  try
    load(fullfile(savdir,'L4alig.mat'));
    L4tot(monkeyIt)=L4(1);
  catch
    clear L4
  end
  
  penNo=str2double(pendir{1}(end-2:end));
  if any(alignment(:,1)==penNo)
    aligOff(monkeyIt,:)=alignment(alignment(:,1)==penNo,2:end);
  else
    aligOff(monkeyIt,:)=0;
  end
  
  numelec=numel(chansel{monkeyIt});
  stats=cell(numelec,1);
  numShanks=ceil(numelec/16);
  shapeTim{monkeyIt}=cell(16,numShanks);
  zCell=cell(numelec,1);
  zCellAlign=zCell;
  for elec=1:numelec
    fname=fullfile(savdir, [fnameBase num2str(chansel{monkeyIt}(elec),'%02d')]);
    loaded=0;
    itcount=1;
    while ~loaded && itcount<1000
      try
        load(fname)
        loaded=1;
      catch
        itcount=itcount+1;
      end
    end
    
    if ~loaded
      error(['loading of ' fname ' failed'])
    end
    
    shapeTim{monkeyIt}{elec}=cfg.best_loc;
    
    fs=cfg.fs;
    cfg.newLen=round(1.5*cfg.winLen);
    [z,s]=bg_swm_extract(cfg);
    
    zsel=~any(isnan(z),2);
    z=z(zsel,:);
    s=s(zsel,:);
    %     c=sum(bsxfun(@minus,z,nanmean(z)).^2,2);
    %     zsel=c<median(c(~isnan(c)));
    
    zCell{elec}=z;
    
    s=bsxfun(@minus,s,nanmean(s,2));
    
    stats{elec}=bg_bootstrap_interpolate(s(:,:),100,1,1,4);
    
    load(cfg.fname, cfg.varname)
    eval(['dat = ' cfg.varname ';'])
    clear(cfg.varname)
    if elec==1
      datCorr=nan(sum(~isnan(dat(:))),numelec);
    end
    nandum=isnan(dat);
    dat(nandum)=0;
    dat=ft_preproc_lowpassfilter(dat,fs,60);
    dat(nandum)=nan;
    dat=dat';
    datCorr(:,elec)=dat(~isnan(dat));
  end
  st=bg_concat_struct(stats);
  skew{monkeyIt}=nan(numel(st),3);
  period=skew{monkeyIt};
  for k=1:numel(st)
    skew{monkeyIt}(k,1)=st(k).skw.mu;
    skew{monkeyIt}(k,2)=st(k).skw.sem;
    skew{monkeyIt}(k,3)=st(k).skw.p_t<(.05/numelec);
    period(k,1)=st(k).period.mu;
    period(k,2)=st(k).period.sem;
    period(k,3)=0;%st(k).period.p_t<(.05/numelec);
  end
  
  [r{monkeyIt},p]=corr(datCorr,'rows','complete');
  
  %% make/align mean shapes
  
  st=reshape(st,[],numShanks);
  shapLen=size(st(1).meanShape,1)*2;
  Shape{monkeyIt}=nan(shapLen,2,numShanks);
  t{monkeyIt}=([1:shapLen]-floor(shapLen/2))/100/fs*1e3;
  mExtr{monkeyIt}=nan(2,2,numShanks);
  for k=1:numShanks
    extr=cat(1,st(:,k).extrema);
    maxsz=max(diff(extr(:,[1 3]),1,2));
    %       [~,pOffSet]=min(abs(bsxfun(@minus,extr(:,2),median(extr))),[],2);
    
    %       offSet=median(extr(:,2))-extr(sub2ind(size(extr),[1:16]',pOffSet));
    %       pOffSet=pOffSet-median(1:size(extr,2));
    [~,offSet]=min([st(:,k).meanShape]);
    
    % align the cfg.best_loc's for cross-correlations
    offSetLoc=offSet/100;
    offSetLoc=offSetLoc-min(offSetLoc);
    winLen=size(zCell{1},2);
    for elec=1:16
      shapeTim{monkeyIt}{elec,k}=shapeTim{monkeyIt}{elec,k}+offSetLoc(elec);
      zCellAlign{elec}=nan(size(zCell{elec})+[0 ceil(max(offSetLoc))]);
      zCellAlign{elec}(:,-round(offSetLoc(elec))+[1:winLen]+ceil(max(offSetLoc)))=zCell{elec};
    end
    
    % do superficial vs deep stats
    zSup=cat(1,zCellAlign{[max(aligOff(k)+1,1):floor(L4tot(monkeyIt))]});
    zDeep=cat(1,zCellAlign{max(ceil(L4tot(monkeyIt)),1):end});
    nanDum=~any(isnan(zSup),1);
    zSup=zSup(:,nanDum);
    zDeep=zDeep(:,nanDum);
    
    mainStats{1,monkeyIt}=bg_bootstrap_interpolate(zSup,100,1,1,3);
    mainStats{2,monkeyIt}=bg_bootstrap_interpolate(zDeep,100,1,1,3);
    sem(:,monkeyIt)=[mainStats{1,monkeyIt}.skw.sem; mainStats{2,monkeyIt}.skw.sem];
    semxy=sqrt(sum(sem(:,monkeyIt).^2));
    mu(:,monkeyIt)=[mainStats{1,monkeyIt}.skw.mu; mainStats{2,monkeyIt}.skw.mu];
    sampsz=[mainStats{1,monkeyIt}.sampSz mainStats{2,monkeyIt}.sampSz];
    t_main=diff(mu(:,monkeyIt))/semxy;
    df=semxy^4/(sem(1)^4/(sampsz(1)-1)+sem(2)^4/(sampsz(2)-1));
    ttest(monkeyIt)=tcdf(-abs(t_main),df);
    
    
    
    aligShapes=nan(shapLen,16);
    for n=1:16
      
      selVec=round([1:size(st(n,k).meanShape,1)]+shapLen/2-offSet(n));
      aligShapes(selVec,n)=st(n,k).meanShape;
    end
    Shape{monkeyIt}(:,1,k)=nanmean(aligShapes(:,[max(aligOff(k)+1,1):floor(L4tot(monkeyIt))]),2);
    Shape{monkeyIt}(:,2,k)=nanmean(aligShapes(:,max(ceil(L4tot(monkeyIt)),1):end),2);
    [~,mExtrdum]=findpeaks(Shape{monkeyIt}(:,1,k));
    [~,selVec]=sort(abs(mExtrdum-size(Shape{monkeyIt},1)/2));
    mExtr{monkeyIt}(:,1,k)=mExtrdum(selVec([1 2]));
    [~,mExtrdum]=findpeaks(Shape{monkeyIt}(:,2,k));
    [~,selVec]=sort(abs(mExtrdum-size(Shape{monkeyIt},1)/2));
    mExtr{monkeyIt}(:,2,k)=mExtrdum(selVec([1 2]));
    

    
%     
  end
  
  %% calculate xcorrs
  % filter kernel
  sig=1;
  h_filt=exp(-([-ceil(5*sig):ceil(5*sig)].^2)./(2*sig).^2);
  if sig==0;
    h_filt=1;
  end
  maxLags=150;
  lags=-maxLags:maxLags;
  numPairs=16*15/2+16;
  pairs=nan(numPairs,2);
  cnt=0;
  xCorrTimeSeries=zeros(2*maxLags+1,numPairs);
  
  xCorrMat=zeros(16,16,numShanks);
  xCorrMatLag=zeros(16,16,numShanks);
  
  for n=1:16
    for m=1:n
      cnt=cnt+1;
      pairs(cnt,:)=[n m];
    end
  end
  
  
  for k=1:numShanks
    numTrials=size(shapeTim{monkeyIt}{1,k},1);
    for n=1:numTrials
      xCorrDum=zeros(ceil(max(shapeTim{monkeyIt}{elec,k}(n,:))),16);
      for elec=1:16
        shapePosIdx=shapeTim{monkeyIt}{elec,k}(n,:);
        shapePosIdx=shapePosIdx(~isnan(shapePosIdx));
        xCorrDum(round(shapePosIdx),elec)=1;
      end
      
      % smooth it
      xCorrDum=filter(h_filt,sum(h_filt),xCorrDum);
      
%       [c(:,n)]=xcorr(mean(xCorrDum(:,9:14),2)-mean(mean(xCorrDum(:,9:14)),2),mean(xCorrDum(:,3:8),2)-mean(mean(xCorrDum(:,3:8),2)),maxLags);
      
      % normalize it
      xCorrNorm=sum(xCorrDum.^2);
           
      for pairIdx=1:numPairs
        dum=xcorr(xCorrDum(:,pairs(pairIdx,1)),xCorrDum(:,pairs(pairIdx,2)),maxLags);
        dum=dum/min(xCorrNorm(pairs(pairIdx,:)));
        xCorrTimeSeries(:,pairIdx)=xCorrTimeSeries(:,pairIdx)+dum;
      end
      
    end
    
    for pairIdx=1:numPairs
      if pairs(pairIdx,1) ~= pairs(pairIdx,2)
      [mVal,mPos]=max(xCorrTimeSeries(:,pairIdx));
      xCorrMat(pairs(pairIdx,1),pairs(pairIdx,2),k)=mVal;
      xCorrMat(pairs(pairIdx,2),pairs(pairIdx,1),k)=mVal;
      xCorrMatLag(pairs(pairIdx,1),pairs(pairIdx,2),k)=lags(mPos);
      xCorrMatLag(pairs(pairIdx,2),pairs(pairIdx,1),k)=-lags(mPos);
      end
    end
    xCorrMat=xCorrMat/numTrials;
    xCorrS{monkeyIt}=xCorrMat;
    xCorrL{monkeyIt}=xCorrMatLag;
  end
  
  
end
  

%%
sh=.05;
mb=.08;
mr=.05;
ml=.08;
sv=.12;
fsz=14;
labsz=24;
brdlabs={'L1','L5'};

col=jet(22);
col=bsxfun(@rdivide,col,sqrt(sum(col.^2,2)));





figure(34)
set(gcf,'position',[100 100 1300 600])
clf

for monkeyIt=1:numel(basedir)
  figure(34)
  shankIdx=ceil(chansel{monkeyIt}((k-1)*16+1)/16);
  shankSel=[1:16]+(k-1)*16;
  
  subaxis(3,10,1+(monkeyIt-1)*5,1,2,2,'sh',sh,'mb',mb,'mr',mr,'mt',.1,'ml',ml,'sv',sv)
  bg_figureLabel(char(64+1+(monkeyIt-1)*4),gca,'fontsize',labsz)
  set(gca,'fontsize',fsz-2)
  aligOffDum=aligOff(shankIdx);
  if isnan(aligOffDum)
    aligOffDum=0;
  end
  x=([1:16]-aligOffDum)*150;
  barhplotsemstar((skew{monkeyIt}(shankSel,1)),skew{monkeyIt}(shankSel,2),skew{monkeyIt}(shankSel,3),col([1:16]+5-aligOffDum,:),x);
  set(gca,'ydir','reverse')
  ylabel('cortical depth (\mum)','fontsize',fsz)
  xlabel('Skweness Index','fontsize',fsz)
  xlim((max(abs(skew{monkeyIt}(:,1)))+(max(skew{monkeyIt}(:,2)))*10)*[-1 1])
%   ylim(deepLims)
  if isnan(aligOff(shankIdx))
    hline(.5*150,'--k','unknown')
  else
    hline(.5*150,'--k',brdlabs{1})
  end
  
  if exist('L4tot','var')
    hline((L4tot(monkeyIt)-aligOffDum)*150,'--k',brdlabs{2})
  end
  
  
  
  subaxis(3,10,2.5+(monkeyIt-1)*5,1,3.5,2)
  %       rdum=nan(16+diff(minmax(aligOff)));
  %       rdum([1:16]+max(aligOff)-aligOff(k),[1:16]+max(aligOff)-aligOff(k))=r(shankSel,shankSel);
  set(gca,'fontsize',fsz-2)
%   imagesc(x,x,xCorrS{monkeyIt})
  hh=imagesc(x,x,r{monkeyIt}(shankSel,shankSel),[floor(min(min(r{monkeyIt}(shankSel,shankSel)))) 1]);
  set(gca,'yticklabel',[])
  xlabel('cortical depth (\mum)','fontsize',fsz)
%   ylabel('cortical depth (\mum)','fontsize',fsz)
  ax=get(gca,'position');
  h=colorbar;
  set(get(h,'title'),'string','Pearson''s r','fontsize',fsz)
%   set(gca,'position',ax)
  bg_figureLabel(char(64+2+(monkeyIt-1)*4),gca,'fontsize',labsz)
%   xlim(deepLims)
%   ylim(deepLims)
  
  if isnan(aligOff(shankIdx))
    vline(150*(.5),'-k','unknown')
    vline(150*(.5),'--w')
    hline(150*(.5),'-k','unknown')
    hline(150*(.5),'--w')
  else
    vline(150*(.5),'-k',brdlabs{1})
    vline(150*(.5),'--w')
    hline(150*(.5),'-k',brdlabs{1})
    hline(150*(.5),'--w')
  end
  
  if exist('L4tot','var')
    hline((L4tot(monkeyIt)-aligOffDum)*150,'-k',brdlabs{2})
    hline((L4tot(monkeyIt)-aligOffDum)*150,'--w')
    vline((L4tot(monkeyIt)-aligOffDum)*150,'-k',brdlabs{2})
    vline((L4tot(monkeyIt)-aligOffDum)*150,'--w')
  end
  subaxis(3,8,(monkeyIt-1)*4+1,3,3,1)
  h=plot(t{monkeyIt},Shape{monkeyIt}(:,:,k));
  set(h,{'color'},{[0 0 .9];[.9 0 0]})
  legend('Superficial','deep')
  xlim([-40 40])
  ylim(maxabs(Shape{monkeyIt})*1.1)
  vline(t{monkeyIt}(mExtr{monkeyIt}(:,1,k)),'--b')
  vline(t{monkeyIt}(mExtr{monkeyIt}(:,2,k)),'--r')
  vline(0,'--k')
  xlabel('time (ms)','fontsize',fsz)
  ylabel('LFP (mV)','fontsize',fsz)
  xlim([-1 1]*30);
  bg_figureLabel(char(64+3+(monkeyIt-1)*4),gca,'fontsize',labsz)
  
  
  subaxis(3,8,(monkeyIt-1)*4+4,3,1,1)
  barplotsemstar(mu(:,monkeyIt),sem(:,monkeyIt),[],[0 0 .9; .9 0 0]);
  ax=axis;  
  dyax=diff(ax([3 4]));
  ax(4)=ax(4)+dyax/10;
  if ttest(monkeyIt)<.001
    xbrack=[1 1 1.5 1.5 1.5 2 2];
    ybrack=[1.5 1.75 1.75 2 1.75 1.75 1.5]*max([mu(:,monkeyIt)+max(sem(:,monkeyIt)); max(sem(:,monkeyIt))]);
    hold on
    plot(xbrack,ybrack,'k')
    axis tight
    ax=axis;
    dyax=diff(ax([3 4]));
    text(1.5,ax(4)+dyax/10,'\ast','fontsize',14,'HorizontalAlignment','center');
    ax(4)=ax(4)+dyax/5;
   
  end
  ax([1 2])=[.5 2.5];
  ax(3)=ax(3)-dyax/10;
  axis(ax)
  ylabel('Skewness Index','fontsize',fsz)
  set(gca,'xtick',[1 2],'xticklabel',{'Sup','Deep'})
  set(gca,'fontsize',fsz-2)
  bg_figureLabel(char(64+4+(monkeyIt-1)*4),gca,'fontsize',labsz)
end

tOff=[-.25 .25];
for monkeyIt=1:2
  ht=mtit(Monkey{monkeyIt},'yoff',-.08,'fontsize',18);
  set(ht.ah,'position',[0+mr 0 1-ml 1])
  tPos=get(ht.th,'position');
  tPos(1)=tPos(1)+tOff(monkeyIt);
  set(ht.th,'position',tPos);
end


%%
export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs/Paper','Fig5_monkeyStats'),'-transparent','-png','-eps','-pdf',34)