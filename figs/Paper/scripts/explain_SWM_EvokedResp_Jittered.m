fs=1e3;
dt=1/fs;
t=0:dt:1;

tt=[0:dt:.5]';

tau=[.2 .1]/4;
f=15;
powerFac=-.5;

signal=(exp(-(tt/tau(1)))-exp(-(tt/tau(2)))).*sin(tt.*((tt).^powerFac)*2*pi*f);
signal(isnan(signal))=0;
signal=signal./(sqrt(sum(signal.^2)));
snr=.7;
%
numTrials=200;
dat=randn(numel(t),numTrials);

[s, f]=bg_mtmfft(dat,fs);
% 1/f filter
kernel=f.^-.6;
kernel(f==0)=0;
% kernel(f<=2)=0;
kernel=kernel;
s=bsxfun(@times,s,kernel);
% s=detrend(s);

dat=real(ifft(s));
dat=ft_preproc_highpassfilter(dat',fs,[2])';
dat=detrend(dat);
dat=zscore(dat)/snr;
signal=zscore(signal);


startPos=randperm(numel(t)-numel(tt),numTrials);

realMean=zeros(size(tt));
for n=1:numTrials
  dat(startPos(n)-1+[1:numel(tt)],n)=dat(startPos(n)-1+[1:numel(tt)],n)+signal;
  realMean=realMean+dat(startPos(n)-1+[1:numel(tt)],n);
end
realMean=realMean/numTrials;

%%
winLen=200;
figure(1)
clf
subaxis(2,1,1)
plot(tt,signal,tt,realMean)
legend('real signal','aligned avg')
xlim([1 winLen]/fs)
subaxis(2,1,2)
plot(t,dat(:,1))
vline(t(startPos(1)))




%%

cfg=[];
cfg.winLen=winLen;
cfg.Tfac=logspace(-4,-1,20);
cfg.numWin=1;
cfg.numIt=2e4;
cfg.fullOutput=1;
cfg.konstant=100;

cfg=bg_SWM(cfg,dat');

%% align using correlation

[r, lags]=bg_windowCorr(cfg.best_s,dat');
[fitVal,midx]=max(r,[],2);

%%
fsz=14;

figure(2)
clf
set(gcf,'position',[100 100 1200 600])
numPlots=4;
trialNo=randperm(numTrials,numPlots);
for n=1:numPlots
  if n==numPlots
    subaxis(numPlots+1,3,1,n+1,2,1)
  else
    subaxis(numPlots+1,3,1,n,2,1,'sh',.1,'mt',.05)
  end
  set(gca,'fontsize',fsz)
  plot(dat(:,trialNo(n)),'k');
  x=[0 winLen-1 winLen-1 0 0]+startPos(trialNo(n));
  y=[-5 -5 5 5 -5];
  ylim([-5 5])
  hold on
  plot(x,y,'-','color',[1 1 1]*.4,'linewidth',1)
  set(gca,'visible','off')
  xlim([1 numel(t)])
  if n==numPlots
    ylabel(['Trial N'])
  else
    ylabel(['Trial ' num2str(n)])
  end
  set(get(gca,'yLabel'),'visible','on')
end


subaxis(numPlots+1,3,1,n,2,1)
set(gca,'fontsize',fsz)
scatter([-1 0 1],[0 0 0],'ok','filled')
xlim([-1 1]*15)
axis off
ylabel('...')
set(get(gca,'yLabel'),'visible','on')

subaxis(3,3,3,2)
set(gca,'fontsize',fsz)
plot(realMean,'k','linewidth',2)
xlim([1 winLen])
x=[1 winLen winLen 1 1];
y=[-5 -5 5 5 -5];
hold on
plot(x,y,'-k','color',[1 1 1]*.4,'linewidth',1)
title('Mean Shape')
axis off

%%
figdir='/home/mrphys/bargip/GIT/SWM/figs/Paper';

export_fig(fullfile(figdir,'explainSWM2'),'-png','-eps','-pdf','-transparent','-nocrop',2)