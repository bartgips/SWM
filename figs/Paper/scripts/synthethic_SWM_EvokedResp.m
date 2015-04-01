clear all
tic
submitBatch=0;
doAnalysis=1;
savdir='/home/mrphys/bargip/GIT/SWM/figs/stats/synthetic_EvokedResp/pinkNoise_fhp5/';
if submitBatch
  
  
  %% loading/compiling function
  
  load ~/GIT/SWM/compiled/latest
  revNumImp=importdata('/home/mrphys/bargip/GIT/SWM/version');
  revNumSC=regexp(revNumImp,'(r\S+-\d+)-','tokens');
  if isempty(revNumSC{1})
    revNumSC=regexp(revNumImp,'(r\S+)','tokens');
  end
  revNumSC=strrep(revNumSC{1}{1}{1},'-','_');
  if ~strcmp(revNumSC,revNum)
    revNum=revNumSC;
    temp=pwd;
    cd ~/GIT/SWM/compiled
    batchid=['bg_SWM_' revNum];
    compiledfun=qsubcompile('bg_SWM','batchid',batchid);
    compiledfun.revNum=revNum;
    fName=fullfile('~/GIT/SWM/compiled',batchid);
    save(fName,'compiledfun','revNum');
    save('~/GIT/SWM/compiled/latest','compiledfun','revNum')
    cd(temp)
  end
  disp(['loaded SWM function version ''' revNum ''''])
  %%
  SNR=[.1:.02:.6];
  
  
  if ~exist(savdir,'file')
    mkdir(savdir)
  end
  
  
  numTrials=200;
  fs=1e3;
  dt=1/fs;
  t=0:dt:1;
  tLen=numel(t);
  winLen=200;
  tt=[0:dt:.5]';
  
  tau=[.3 .1]/4;
  baseFreq=10;
  powerFac=.6;
  
  signal=(exp(-(tt/tau(1)))-exp(-(tt/tau(2)))).*sin(tt.^powerFac.*2*pi*baseFreq);
  signal(isnan(signal))=0;
  signal=signal./(sqrt(sum(signal.^2)));
  signal=zscore(signal);
  
  if exist(fullfile(savdir,'syntDat.mat'),'file')
    load(fullfile(savdir,'syntDat.mat'),'startPos')
    save(fullfile(savdir,'syntDat.mat'),'SNR','-append');
  else
    startPos=randperm(numel(t)-numel(tt),numTrials);
    save(fullfile(savdir,'syntDat.mat'),'SNR','startPos','signal');
  end
  
  
  
  
  if ~exist(fullfile(savdir,'torque'),'file')
    mkdir(fullfile(savdir,'torque'))
  end
  cd(fullfile(savdir,'torque'))
  for snrIdx=1:numel(SNR)
    
    snr=SNR(snrIdx);
    %
    
    dat=randn(numel(t),numTrials);
    
    [s, f]=bg_mtmfft(dat,fs);
    % 1/f filter
    kernel=abs(f).^-1+.01;
    kernel(f==0)=0;
    %   kernel(abs(f)<baseFreq/2)=0;
    kernel(f<=1)=0;
    kernel=kernel;
    s=bsxfun(@times,s,kernel);
    %   s=detrend(s);
    
    dat=real(ifft(s));
    % dat=ft_preproc_highpassfilter(dat',fs,[2])';
    % dat=detrend(dat);
    dat=zscore(dat)/snr;
    [ssig,fsig]=bg_mtmfft(signal,fs,0,0,0,size(dat,1));
    [s,f]=bg_mtmfft(dat,fs);
    
    
    
    
    
    
    
    realMean=zeros(size(tt));
    for n=1:numTrials
      dat(startPos(n)-1+[1:numel(tt)],n)=dat(startPos(n)-1+[1:numel(tt)],n)+signal;
      realMean=realMean+dat(startPos(n)-1+[1:numel(tt)],n);
    end
    realMean=realMean/numTrials;
    
    [s2,f2]=bg_mtmfft(dat,fs);
    %
    %   figure(111)
    %   plot(f(f>0),mean(mean(abs(s(f>0,:,:)).^2,3),2),fsig(fsig>0),abs(ssig(fsig>0)).^2,f2(f2>0),mean(mean(abs(s2(f2>0,:,:)).^2,3),2));
    %   legend('noise','signal','sum')
    %   xlim([0 100])
    %
    %
    %
    %   figure(1)
    %   clf
    %   subaxis(2,1,1)
    %   plot(tt,signal,tt,realMean)
    %   legend('real signal','aligned avg')
    %   xlim([1 winLen]/fs)
    %   subaxis(2,1,2)
    %   idx=randperm(numTrials,1);
    %   plot(t,dat(:,idx))
    %   vline(t(startPos(idx)))
    
    
    
    
    
    %%
    
    
    cfg=[];
    cfg.winLen=winLen;
    cfg.Tfac=logspace(-4,-1,30);
    cfg.winPerTrial=1;
    cfg.numIt=5e5;
    cfg.fullOutput=0;
    cfg.konstant=100;
    cfg.verbose=0;
    cfg.Fhp=8;
    cfg.fs=fs;
    cfg.revNum=revNum;
    % cfg.best_loc=startPos(:);
    
    fName=['cfg_' num2str(snrIdx,'%03d')];
    
    varName=['dat' num2str(snrIdx,'%03d')];
    
    eval([varName ' = dat'';'])
    
    
    
    cfg.fName=fullfile(savdir,'syntDat.mat');
    cfg.varName=varName;
    cfg.outputFile=fullfile(savdir, fName);
    
    
    batchid=fName;
    
    
    if ~exist(cfg.outputFile,'file')  && ~exist([cfg.outputFile '.mat'],'file')
      save(fullfile(savdir,'syntDat.mat'),varName,'-append')
      qsubfeval(compiledfun, cfg, 'memreq', 1*1024^3, 'timreq', 2*3600, 'batchid', [batchid]);
    else
      warning(sprintf(['File ' cfg.outputFile ' already exists.\nNot re-running SWM...']))
    end
    
    
  end
end


if doAnalysis
  %% Analysis
  
  outpDir=savdir;
  load(fullfile(outpDir,'syntDat.mat'),'SNR','startPos','signal');
  
  % SNR=SNR(1:17);
  cmap=jet(numel(SNR));
  
  snrlims=[min(SNR) min(1,max(SNR))]+min(diff(SNR))*.5*[-1 1];
  
  mShapesS=zeros(250,numel(SNR));
  mShapesZ=mShapesS;
  histBins=-1000:1000;
  N=nan(numel(histBins),numel(SNR));
  for snrIdx=1:numel(SNR)
    isLoaded=0;
    while ~isLoaded
      try
        load(fullfile(outpDir,['cfg_' num2str(snrIdx,'%03d')]))
        isLoaded=1;
      end
    end
    mShapesS(1:cfg.winLen,snrIdx)=cfg.best_s;
    mShapesZ(1:cfg.winLen,snrIdx)=cfg.best_z;
    N(:,snrIdx)=hist(cfg.best_loc(:)-startPos(:),histBins);
  end
  
  %% align all found shapes
  
  for snrIdx=1:numel(SNR)
    load(fullfile(outpDir,['cfg_' num2str(snrIdx,'%03d')]))
    [xc,lags]=xcorr(mShapesZ(:,snrIdx),signal(:));
    [~,midx]=max(xc);
    cfg.best_loc=cfg.best_loc+lags(midx);
    cfg.winLen=cfg.winLen+50;
    cfg=rmfield(cfg,'Fhp');
    [S,Z]=bg_swm_extract(cfg);
    mShapesZ(:,snrIdx)=nanmean(Z);
    mShapesS(:,snrIdx)=nanmean(S);
    N(:,snrIdx)=hist(cfg.best_loc(:)-startPos(:),histBins);
  end
  
%% use xCorr with data to find shapes based on mean shape
  
  meanDist=nan(numel(SNR),1);
  stdDist=meanDist;
  varDist=meanDist;
  
  for snrIdx=1:numel(SNR)
    load(fullfile(outpDir,['cfg_' num2str(snrIdx,'%03d')]))
    cfg.winLen=cfg.winLen+50;
    datname=num2str(snrIdx,'dat%03d');
    load(fullfile(outpDir,'syntDat'),datname)
    eval(['dat=' datname ';'])
    dat=ft_preproc_highpassfilter(dat,cfg.fs,cfg.Fhp);
    clear(datname)
    [r,lags]=bg_windowCorr(mShapesZ(:,snrIdx),dat);
    [~,midx]=max(r,[],2);
    cfg.best_loc=lags(midx)';
    [S,Z]=bg_swm_extract(cfg,dat);
    % one more alignment
    [xc,lags]=xcorr(nanmean(Z),signal);
    [~,midx]=max(xc);
    cfg.best_loc=cfg.best_loc+lags(midx);
    cfg=rmfield(cfg,'Fhp');
    [S,Z]=bg_swm_extract(cfg,dat);
    mShapesZ(:,snrIdx)=nanmean(Z);
    mShapesS(:,snrIdx)=nanmean(S);
    N(:,snrIdx)=hist(cfg.best_loc(:)-startPos(:),histBins);
    meanDist(snrIdx)=mean(abs(cfg.best_loc(:)-startPos(:)));
    stdDist(snrIdx)=std(abs(cfg.best_loc(:)-startPos(:)));
    distrib=cfg.best_loc(:)-startPos(:);
    [bs, distr]=bg_bootstrap_std(distrib,1,1e3);
%     varDist(snrIdx)=var(cfg.best_loc(:)-startPos(:));
    varDist(snrIdx)=bs.mu^2;
%     CIDist(snrIdx,:)=bs.CI;
    [bs, distr]=bg_bootstrap_mean(abs(distrib),1,1e3);
    meanDist(snrIdx)=bs.mu;
  end
  
  
  
  
  %% ffts
  fs=1e3;
  load(fullfile(outpDir,'syntDat'),'signal')
  nfft=4*fs+1;
  datFT=nan(nfft,numel(SNR));
  
  for snrIdx=1:numel(SNR)
    datname=num2str(snrIdx,'dat%03d');
    load(fullfile(outpDir,'syntDat'),datname)
    eval(['dat=' datname ';'])
    clear(datname)
    [ftdum,f]=bg_mtmfft(dat.',fs,0,0,0,nfft);
    datFT(:,snrIdx)=nanmean(nanmean(abs(ftdum).^2,3),2);
  end
  
  
  %% plot results
  %interpolate (if necessary)
  SNRi=[min(SNR):min(diff(SNR)):max(SNR)];
  mZi=interp1(SNR(:),mShapesZ.',SNRi(:));
  mSi=interp1(SNR(:),mShapesS.',SNRi(:));
  FTi=interp1(SNR(:),datFT.',SNRi(:));
  Ni=interp1(SNR(:),N.',SNRi(:));
  
  freqlims=[0 40];
  
  figure(10)
  set(gcf,'position',[100 100 1200 800])
  clf
  subaxis(3,2,1,1,'sv',.08)
  imagesc(f(f>=0),SNRi,bsxfun(@rdivide,FTi(:,f>=0),mean(FTi,2)))
  colorbar
  xlabel('Freq (Hz)')
  ylabel('SNR')
  title('normalized PSD')
  xlim(freqlims)
  subaxis(3,2,2,1,'sv',.08)
  imagesc(f(f>=0),SNRi,FTi(:,f>=0))
  colorbar
  xlabel('Freq (Hz)')
  ylabel('SNR')
  ylim(snrlims)
  title('raw PSD')
  % imagesc(f(f>=0),SNR,datFT(f>=0,:)')
  xlim(freqlims)
  subaxis(3,2,1,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mZi)
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('found mean shape (zscored)')
  subaxis(3,2,2,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mSi)
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('found mean shape (raw)')
  
  subaxis(3,1,3)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  imagesc(histBins,SNRi,Ni/cfg.numWindows)
  xlim([-20 20])
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('devation from true position of evoked response')
  title('probability histogram window location')
  
  %% fourier equivalent
  fs=1e3;
  load(fullfile(outpDir,'syntDat'),'signal')
  signal=zscore(signal);
  [signalfft,ff]=bg_mtmfft(signal(:),fs,0,0,0,nfft);
  signalfft=signalfft(ff>=0,:);
  ff=ff(ff>=0);
  %find bandpass
  signalPow=abs(signalfft).^2;
  maxPow=max(signalPow);
  powFac=.1;
  fbp=find(signalPow>powFac*maxPow,1,'first');
  fbp(2)=find(signalPow>powFac*maxPow,1,'last');
  fbp=ff(fbp);
  
  figure(13)
  plot(ff,abs(signalfft).^2)
  vline(fbp)
  xlim(freqlims)
  
  mShapesSft=nan(250,numel(SNR));
  mShapesZft=mShapesS;
  histBins=-1000:1000;
  Nft=nan(numel(histBins),numel(SNR));
  meanDistft=nan(numel(SNR),1);
  stdDistft=meanDistft;
  varDistft=meanDistft;
  
  for snrIdx=1:numel(SNR)
    datname=num2str(snrIdx,'dat%03d');
    load(fullfile(outpDir,'syntDat'),datname)
    eval(['dat=' datname ';'])
    clear(datname)
    datFilt=ft_preproc_bandpassfilter(dat,fs,fbp);
    ampEnv=abs(hilbert(datFilt.'));
    [~,midx]=max(ampEnv);
    load(fullfile(outpDir,['cfg_' num2str(snrIdx,'%03d')]))
    cfg.best_loc=midx(:);
    cfg.winLen=size(mShapesZft,1);
    [S,Z]=bg_swm_extract(cfg,dat);
    [r,lags]=bg_windowCorr(nanmean(Z)',dat);
    [~,midx]=max(r,[],2);
    cfg.best_loc=lags(midx)';
    [S,Z]=bg_swm_extract(cfg,dat);
    % one more alignment
    [xc,lags]=xcorr(nanmean(Z),signal);
    [~,midx]=max(xc);
    cfg.best_loc=cfg.best_loc+lags(midx);
    cfg=rmfield(cfg,'Fhp');
    [S,Z]=bg_swm_extract(cfg,dat);
    mShapesZft(:,snrIdx)=nanmean(Z);
    mShapesSft(:,snrIdx)=nanmean(S);
    Nft(:,snrIdx)=hist(cfg.best_loc(:)-startPos(:),histBins);
    meanDistft(snrIdx)=mean(abs(cfg.best_loc(:)-startPos(:)));
    stdDistft(snrIdx)=std(abs(cfg.best_loc(:)-startPos(:)));
    varDistft(snrIdx)=var(cfg.best_loc(:)-startPos(:));
    distrib=cfg.best_loc(:)-startPos(:);
    [bs, distr]=bg_bootstrap_std(distrib,1,1e3);
    varDistft(snrIdx)=bs.mu^2;
    [bs, distr]=bg_bootstrap_mean(abs(distrib),1,1e3);
    meanDistft(snrIdx)=bs.mu;
  end
  
  %% plot ft results
  %interpolate (if necessary)
  SNRi=[min(SNR):min(diff(SNR)):max(SNR)];
  mZift=interp1(SNR(:),mShapesZft.',SNRi(:));
  mSift=interp1(SNR(:),mShapesSft.',SNRi(:));
  FTi=interp1(SNR(:),datFT.',SNRi(:));
  Nift=interp1(SNR(:),Nft.',SNRi(:));
  
  
  
  figure(11)
  set(gcf,'position',[100 100 1200 800])
  clf
  subaxis(3,2,1,1,'sv',.08,'mt',.05,'mr',.05,'ml',.05)
  imagesc(f(f>=0),SNR,bsxfun(@rdivide,FTi(:,f>=0),mean(FTi,2)))
  colorbar
  xlabel('Freq (Hz)')
  ylabel('SNR')
  ylim(snrlims)
  title('normalized PSD')
  xlim(freqlims)
  subaxis(3,2,2,1,'sv',.08)
  imagesc(f(f>=0),SNR,FTi(:,f>=0))
  colorbar
  xlabel('Freq (Hz)')
  ylabel('SNR')
  ylim(snrlims)
  title('raw PSD')
  % imagesc(f(f>=0),SNR,datFT(f>=0,:)')
  xlim(freqlims)
  subaxis(3,2,1,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mZift)
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('found mean shape (zscored)')
  subaxis(3,2,2,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mSift)
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('found mean shape (raw)')
  
  subaxis(3,1,3)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  imagesc(histBins,SNRi,Nift/cfg.numWindows)
  xlim([-20 20])
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('devation from true position of evoked response')
  title('probability histogram window location')
  
  
  %% combined figure
  cmap=jet(2)*2;
  
  figure(12)
  set(gcf,'position',[100 100 1200 800])
  clf
  subaxis(3,2,1,1,'sv',.08,'mt',.05,'mr',.05,'ml',.05)
  set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % imagesc(f(f>=0),SNR,bsxfun(@rdivide,FTi(:,f>=0),max(FTi,[],2)))
  plot(f(f>=0),bsxfun(@rdivide,FTi([1 end],f>=0),max(FTi([1 end],f>=0),[],2)),'linewidth',2)
  % colorbar
  legend(['SNR = ' num2str(SNR(1))],['SNR = ' num2str(SNR(end))])
  xlabel('Freq (Hz)')
  ylabel('Normalized power (a.u.)')
  % ylabel('SNR')
  % ylim(snrlims)
  ylim([0 .5])
  xlim([1 40])
  title('normalized PSD')
  % xlim(freqlims)
  subaxis(6,2,2,1)
  plot((signal),'linewidth',2)
  xlim([0 250])
  title('Evoked response shape')
  xlabel('Time (ms)')
  ylabel('signal (a.u.)')
  subaxis(6,2,2,2)
  % imagesc(f(f>=0),SNR,FTi(:,f>=0))
  % colorbar
  plot(ff,signalPow,'linewidth',2)
  h=vline(fbp,'-k');
  legend(h(1),'filter BW')
  set(h,'linewidth',2)
  xlabel('Freq (Hz)')
  ylabel('Power (a.u.)')
  % ylabel('SNR')
  % ylim(snrlims)
  title('PSD of evoked response')
  % imagesc(f(f>=0),SNR,datFT(f>=0,:)')
  xlim(freqlims)
  subaxis(3,2,1,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mSi,maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]))
  % imagesc(1,SNRi,mZi,maxabs([mZi(SNRi>.3,:) mZift(SNRi>.3,:)]))
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('Found mean shape (SWM)')
  subaxis(3,2,2,2)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  % plot(mShapesZ,'linewidth',2)
  imagesc(1,SNRi,mSift,maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]))
  % imagesc(1,SNRi,mZift,maxabs([mZi(SNRi>.3,:) mZift(SNRi>.3,:)]))
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('Time (ms)')
  title('Found mean shape (PA)')
  
  subaxis(3,3,1,3)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  imagesc(histBins,SNRi,Ni/cfg.numWindows,[0 max([Ni(:); Nift(:)])/cfg.numWindows])
  xlim([-20 20])
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('devation from true position of evoked response (ms)')
  title('probability histogram detected evoked response location (SWM)')
  
  subaxis(3,3,2,3)
  % set(gca,'ColorOrder',cmap,'nextPlot','replacechildren')
  imagesc(histBins,SNRi,Nift/cfg.numWindows,[0 max([Ni(:); Nift(:)])/cfg.numWindows])
  xlim([-20 20])
  colorbar
  ylabel('SNR')
  ylim(snrlims)
  xlabel('devation from true position of evoked response (ms)')
  title('probability histogram detected evoked response location (Fourier)')
  
  subaxis(3,3,3,3)
  cla
  h1=plot(SNR,sqrt(varDist),'b','linewidth',2);
  hold
  h2=plot(SNR,sqrt(varDistft),'r','linewidth',2);
%   h1=plot(SNR,[meanDist meanDist-stdDist/sqrt(cfg.numWindows) meanDist+stdDist/sqrt(cfg.numWindows)],'b');
%   hold on
%   h2=plot(SNR,[meanDistft meanDistft-stdDistft/sqrt(cfg.numWindows) meanDistft+stdDistft/sqrt(cfg.numWindows)],'r');
  legend([h1(1),h2(1)],'SWM','Fourier')
  xlabel('SNR')
  ylabel('Std of error (ms)')
%   ylabel('mean error (ms)')

%% paper figure

fsz=12;

ylims=[-1 1]*5;
snrSel=floor(linspace(1,numel(SNR),6));
snrSel=snrSel(2:end);
cmap=[linspace(0,1,numel(snrSel))' zeros(numel(snrSel),1) linspace(1,0,numel(snrSel))'];

figure(121)
set(gcf,'position',[100 100 600 600])
clf
subaxis(2,2,1,'sv',.08,'sh',.12)
imagesc(1,[SNRi .62 .64],[mSi; nan(1,250); signal(1:250)'],maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]))

im=real2rgb([mSi; nan(1,250); signal(1:250)'],'jet',maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]));
im(end-1,:,:)=1;

colorbar
ylabel('SNR','fontsize',fsz)
% ylim(snrlims)
xlabel('Time (ms)','fontsize',fsz)
title('Mean shape (SWM)','fontsize',fsz)
set(gca,'ytick',[.1:.1:.6 .64],'yticklabel',[num2str([.1:.1:.6]'); 'inf'])
hold on
image(1,[SNRi .62 .64],im)
bg_figureLabel('A')

subaxis(2,2,2)
imagesc(1,[SNRi .62 .64],[mSift; nan(1,250); signal(1:250)'],maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]))
hcb=colorbar;

im=real2rgb([mSift; nan(1,250); signal(1:250)'],'jet',maxabs([mSi(SNRi>.3,:) mSift(SNRi>.3,:)]));
im(end-1,:,:)=1;


% set(get(hcb,'ylabel'),'string','error in found shape')
% ylabel('SNR')
% ylim(snrlims)
xlabel('Time (ms)','fontsize',fsz)
title('Mean shape (Fourier)','fontsize',fsz)
set(gca,'ytick',[.1:.1:.6 .64],'yticklabel',[num2str([.1:.1:.6]'); 'inf'])
pause(.1)
bg_figureLabel('B')
hold on
image(1,[SNRi .62 .64],im)

subaxis(2,2,3)
semilogy(SNR,(sum(bsxfun(@minus,mShapesS,signal(1:250)).^2)),'r','linewidth',2)
hold on
semilogy(SNR,(sum(bsxfun(@minus,mShapesSft,signal(1:250)).^2)),'b','linewidth',2)

bg_figureLabel('C')
legend('SWM','Fourier')
ylabel('Sum of squared error','fontsize',fsz)
xlabel('SNR','fontsize',fsz)

subaxis(2,2,4)
plot(SNR,meanDist,'r','linewidth',2)
hold on
plot(SNR,meanDistft,'b','linewidth',2)
bg_figureLabel('D')
legend('SWM','Fourier')
ylabel(sprintf('mean error in\nwindow location (ms)'),'fontsize',fsz)
xlabel('SNR','fontsize',fsz)

saveFig=0;
if saveFig
figdir='/home/mrphys/bargip/GIT/SWM/figs/Paper';

export_fig(fullfile(figdir,'SimulateERP2'),'-png','-eps','-transparent',121)
end

%% paper figure NEW

fsz=12;

ylims=[-1 1]*5;

labOffset=[-10 -5];

SNRoi=0.3;
SNRSel=bg_findnearest(SNRi,SNRoi);
figure(122)
set(gcf,'position',[100 100 600 600])
clf
subaxis(2,1,1,'sv',.1,'sh',.12)
plot(mSi(SNRSel,:),'r','linewidth',2)
hold on
plot(mSift(SNRSel,:),'linewidth',2)
plot(signal(1:250),'--k','linewidth',2)
xlabel('Time (ms)','fontsize',fsz)
title(['Motifs at SNR=' num2str(SNRoi)],'fontsize',fsz)
bg_figureLabel('A',gca,'offSet',labOffset)
h=legend('SWM','Fourier','original');
bg_shortLegend(h,.5,10);

subaxis(2,2,3)
semilogy(SNR,(sum(bsxfun(@minus,mShapesS,signal(1:250)).^2)),'r','linewidth',2)
hold on
semilogy(SNR,(sum(bsxfun(@minus,mShapesSft,signal(1:250)).^2)),'b','linewidth',2)

bg_figureLabel('B',gca,'offSet',labOffset)
h=legend('SWM','Fourier');
bg_shortLegend(h,.5,10);
ylabel('Sum of squared error','fontsize',fsz)
xlabel('SNR','fontsize',fsz)

subaxis(2,2,4)
plot(SNR,meanDist,'r','linewidth',2)
hold on
plot(SNR,meanDistft,'b','linewidth',2)
bg_figureLabel('C',gca,'offSet',labOffset)
h=legend('SWM','Fourier');
bg_shortLegend(h,.5,10);
ylabel(sprintf('mean error in\nwindow location (ms)'),'fontsize',fsz)
xlabel('SNR','fontsize',fsz)

saveFig=0;
if saveFig
  figdir='/home/mrphys/bargip/GIT/SWM/figs/Paper';
  
  export_fig(fullfile(figdir,'SimulateERP_singlePlot2'),'-png','-eps','-transparent',122)
end
end
toc