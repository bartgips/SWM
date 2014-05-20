
skew=1;
tLen=1000;
verbose=1;


snr=[.1 1 inf];
skw=nan(2,numel(snr));
stat=cell(numel(snr),1);

for n=1:numel(snr)
[skw(:,n),stat{n}]=bg_sawtooth_shapefind_snr(snr(n),skew,tLen,verbose);
end

savedata=0;
if savedata
save fig1data_3 skw stat
end
%%
load fig1data_3
t=stat{1}.t;
signal=[stat{1}.signalOrig stat{2}.signalOrig stat{3}.signalOrig];
signalft=(fft(signal));
nfft=size(signal,1);
signalft=signalft(1:nfft/2+1,:)/nfft;
f=stat{1}.cfgSWM.fs/2*linspace(0,1,nfft/2+1);

sigPow=squeeze(abs(signalft).^2);

meanS_PA=nan(size(stat{1}.meanShape,1),numel(stat));
meanS_SWM=meanS_PA;
for n=1:numel(stat)
  meanS_PA(:,n)=stat{n}.meanShape(:,1);
  meanS_SWM(:,n)=stat{n}.meanShape(:,2);
end
tmean=repmat([1:size(meanS_PA,1)]'/stat{1}.cfgSWM.fs*1e3,1,numel(stat));

%% align mean shapes from SWM

cfg=[];
cfg.winLen=150;
cfg.numWin=1;
cfg.numIt=100;
cfg.Tfac=1e-3;
cfg.verbose=1;
cfg.dispPlot=0;
cfg=bg_SWM_SA(cfg,meanS_SWM');

mShift=cfg.best_loc;
mShift=mShift-mean(mShift);



%%
figure(1)
clf
set(gcf,'position',[50 50 800 600])
subaxis(3,2,1,'sv',.08)

h=plot(t,signal,'linewidth',2);
xlim([0 .5])
set(h,{'color'},num2cell([0 0 1; 1 0 0; 0 0 0],2))
xlabel('time (s)')
ylabel('Signal amplitude')
% title('Simulated Signals')
text(-0.1,1.2,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',24)

subaxis(3,2,2,1)
h=plot(f,sigPow,'linewidth',2);
set(h,{'color'},num2cell([0 0 1; 1 0 0; 0 0 0],2))
xlim([0 25])
legend(num2str(snr'))
xlabel('frequency (Hz)')
ylabel('Power (a.u.)')
ylim([0 max(max(sigPow(f<5,:)))/2])
% title('Power Spectrum')
text(-0.1,1.2,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',24)

subaxis(3,1,2)
h=plot(tmean,fliplr(meanS_PA),'linewidth',2);
% title('Mean shape found by PA')
set(h,{'color'},num2cell(flipud([0 0 1; 1 0 0; 0 0 0]),2))
ylim(maxabs(meanS_PA)*1.1)
text(-0.1,1.2,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',24)
ylabel('Signal amplitude')

subaxis(3,1,3)
h=plot(bsxfun(@minus,tmean,fliplr(mShift')),fliplr(meanS_SWM),'linewidth',2);
% title('Mean shape found by SWM')
set(h,{'color'},num2cell(flipud([0 0 1; 1 0 0; 0 0 0]),2))
xlabel('time (ms)')
xlim(minmax(tmean))
ylim(maxabs(meanS_SWM)*1.1)
ylabel('Signal amplitude')
text(-0.1,1.2,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',24)

export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs','Fig1'),'-transparent','-pdf','-png','-eps');