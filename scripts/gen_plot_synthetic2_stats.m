cd /home/mrphys/bargip/GIT/SWM/figs/stats

%%
skew=1;
tLen=1000;
verbose=0;


snr=[.01:.01:.1 .2:.1:1];
% stat=cell(numel(snr),1);
batchid='Synt_Bat';
numIt=20;

jobid=cell(numIt,1);
for iter=1:numIt
  jobid{iter}=qsubfeval('bg_sawtooth_shapefind_snr_bat',snr,skew,tLen,1,'memreq',1024^3,'timreq',60^2*20,'batchid',batchid);
% [skw(iter,n,:)]=bg_sawtooth_shapefind_snr(snr(n),skew,tLen,verbose);

end

save jobid_1skew_PT_5e3fs jobid snr

%%
load jobid

%%%%%change this

skw=nan([numel(snr) numel(jobid) 2]);
meanShape=nan([numel(snr) numel(jobid) 200 2]);
%%%%%%%
missed=[];
for n=1:numel(jobid)
  home
  disp(['loading result ' num2str([n numel(jobid)],'%d/%d')])
  try
    argout=bg_qsubget(jobid{n});
    skw(:,n,:)=squeeze(argout{1});
    meanShape(:,n,:,:)=squeeze(argout{2});
  catch
    missed=[missed n];
  end
end

save('run1000s_.01_1_1.0.mat', 'skw', 'meanShape', 'missed', 'snr')

%% if satisfied, cleanup
for n=1:numel(jobid)
  qsubget(jobid{n});
end
%% figure
figure(2)
set(gcf,'position',[50 50 1600 800])
clf
h(3)=plot(snr,ones(size(snr))*(25/100-.5)*-2,'--','color',[1 1 1]*.4,'linewidth',2);
hold on
h(1)=plot(snr,nanmean(skw(:,:,1),2),'r','linewidth',2);
plot(snr,quantile(skw(:,:,1),.25,2),'r')
plot(snr,quantile(skw(:,:,1),.75,2),'r')
h(2)=plot(snr,nanmean(skw(:,:,2),2),'linewidth',2);
plot(snr,quantile(skw(:,:,2),.25,2))
plot(snr,quantile(skw(:,:,2),.75,2))
legend(h,'PA','SWM','True Skewness','location','southeast')
xlim(minmax(snr))
xlabel('SNR')
ylabel('Skewness Index')
ylim([0 1])
% export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs','Fig2'),'-transparent','-pdf');