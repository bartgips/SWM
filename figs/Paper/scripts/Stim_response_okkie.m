load('/home/mrphys/bargip/OrientationDataMonkeys/Torque/Okkie/StimResp/cfg_Stim_z_300WL_nW1.mat')
load(cfg.fname,cfg.varname);

%% aligning to a minumum where stimResp is centred in the window
cfg2=cfg;
cfg2.best_loc=cfg.best_loc+110;
[~, z]=bg_swm_extract(cfg2);
zz=bsxfun(@minus,z,nanmean(z));
costOld=sum(nanmean(zz(:,:).^2,2));
zMean=squeeze(nanmean(z));

[r, lags]=bg_windowCorr(zMean, dat_preproc);
[fitVal, midx]=max(r,[],2);
cfg2.best_loc=midx(:);

[~, z]=bg_swm_extract(cfg2);
zz=bsxfun(@minus,z,nanmean(z));
costNew=sum(nanmean(zz(:,:).^2,2));

%%
figure(100)
bar([costOld,costNew])
ylim(minmax([costOld,costNew]).*[.99 1.01])
xlim([0 3])
set(gca,'xtickLabel',{'Old Cost', 'New Cost'})

% export_fig('Stim_z_300WL_nW1_costComp','-png',100)




%%
tSel=[-200:450];
stimLocdum=squeeze(nanmean(dat_preproc(:,500+tSel,:)))';
stimLocdum2=squeeze(nanmean(dat_preproc(:,:,:)))';



cfg2.newLen=400;
[swmLocdum,~,cfg2]=bg_swm_extract(cfg2);
swmLocdum=squeeze(nanmean(swmLocdum))';

[rFit,lags]=bg_windowCorr(swmLocdum',permute(stimLocdum2,[3 2 1]));

[~, fitIdx]=max(rFit);
stimLocdum2=squeeze(nanmean(dat_preproc(:,fitIdx+[1:cfg2.newLen]-1,:)))';


diffIm=swmLocdum-stimLocdum2;

%%
clims=maxabs([maxabs(swmLocdum) maxabs(stimLocdum)]);
fsz=14;
histLims=[-50 30];
contactStart=16;
offSet=[-10 0];

figure(99)
clf
set(gcf,'position',[100 100 1200 400])
subaxis(3,3,1,1,1,2,'sh',.08,'sv',.15)
set(gca,'fontsize',fsz)
imagesc(tSel,150*([1:size(stimLocdum,1)]-contactStart-1),stimLocdum,clims)
xlim([1 cfg2.newLen]-500+fitIdx)
vline(0,'--k')
xlabel('Time since stimulus onset (ms)')
ylim([1 16]*150-150)
ylabel('Depth (\mum)')
h=colorbar;
set(h,'fontsize',fsz)
set(get(h,'ylabel'),'string','LFP (mV)','fontsize',fsz)
title('Aligned to Stimulus onset')
bg_figureLabel('A',gca,'offSet',offSet)

subaxis(3,3,2,1,1,2)
set(gca,'fontsize',fsz)
imagesc([1:cfg2.newLen]-cfg2.addLen,150*([1:size(stimLocdum,1)]-contactStart-1),swmLocdum,clims)
h=colorbar;
set(h,'fontsize',fsz)
set(get(h,'ylabel'),'string','LFP (mV)','fontsize',fsz)
vline([1 300],'-k')
ylim([1 16]*150-150)
xlabel('Time since Sliding Window onset (ms)')
title('SWM')
bg_figureLabel('B',gca,'offSet',offSet)
vline(-cfg2.addLen-fitIdx+500,'--k')

subaxis(3,3,3,1,1,2)
set(gca,'fontsize',fsz)
imagesc([1:cfg2.newLen]-500+fitIdx,150*([1:size(stimLocdum,1)]-contactStart-1),diffIm,maxabs(diffIm))
h=colorbar;
set(h,'fontsize',fsz)
set(get(h,'ylabel'),'string','\DeltaLFP (mV)','fontsize',fsz)
vline(0,'--k')
ylim([1 16]*150-150)
xlabel('Time since stimulus onset (ms)')
title('SWM - Stim Aligned')
bg_figureLabel('C',gca,'offSet',offSet)


subaxis(3,2,1,3,2,1)
set(gca,'fontsize',fsz)
% bar(0,1)
% xlim(histlims)
% subaxis(5,2,2,5,1,1)
[n,x]=hist(cfg2.best_loc-500,[1:1001]-500);
bar(x,n,'k')
xlim(histLims)
xlabel('Window start time relative to stimulus onset (ms)')
ylabel('N')
bg_figureLabel('D',gca,'offSet',offSet)
%%
%%

figure(999)
subaxis(1,2,1)
imagesc(tSel(1:end-1),[1:size(stimLocdum,1)]-contactStart,diffIm,maxabs(diffIm))
colorbar
ylim([1 16])
xlabel('Window start time relative to stimulus onset (ms)')
title('difference between SWM and Stim Locked')
subaxis(1,2,2)
imagesc(tSel(1:end-1),[1:size(stimLocdum,1)]-contactStart,diffIm.*sign(swmLocdum),maxabs(diffIm))
ylim([1 16])
colorbar


export_fig('Stim_z_300WL_nW1_diffLFP','-png',999)

%%
figdir='/home/mrphys/bargip/GIT/SWM/figs/Paper';

export_fig(fullfile(figdir,'Okkie009_StimResp_300ms'),'-png','-eps','-transparent',99)




%% dft filtering 50 Hz
fs=1e3;
t=1:size(dat_preproc,2);
DFTfilt=exp(1i*50*t/fs*2*pi);
% DFTfilt=repmat(DFTfilt,1,48);

for n=1:size(dat_preproc,1)
  datdum=squeeze(dat_preproc(n,:,:));
  fit=datdum'/DFTfilt;
  datdum=(real(datdum'-fit*DFTfilt));
  dat_preproc(n,:,:)=datdum';
end

[s,z]=bg_swm_extract(cfg2);
orig=squeeze(nanmean(z));
[s,z]=bg_swm_extract(cfg2,dat_preproc);
new=squeeze(nanmean(z));

