cd /home/mrphys/bargip/GIT/SWM/figs/stats

%%
skew=1;
tLen=1000;
verbose=0;


snr=[.01:.01:.1 .2:.1:1];
stat=cell(numel(snr),1);

numIt=50;

jobid=cell(numIt,numel(snr));
for iter=1:numIt
for n=1:numel(snr)
  jobid{iter,n}=qsubfeval('bg_sawtooth_shapefind_snr',snr(n),skew,tLen,verbose,'memreq',1024^2*100,'timreq',60*20);
% [skw(iter,n,:)]=bg_sawtooth_shapefind_snr(snr(n),skew,tLen,verbose);
end
end

save jobid jobid snr

%%
load jobid

stats=cell(size(jobid));

skw=nan([size(jobid) 2]);
missed=[];
for n=1:numel(jobid)
  try
  [a,b]=ind2sub(size(jobid),n);
  argout=bg_qsubget(jobid{n});
  skw(a,b,:)=argout{1};
  stats{n}=argout{2};
  catch
    missed=[missed b];
  end
end

save('run100s_.1_1.mat','skw','missed')

%% if satisfied, cleanup
for n=1:numel(jobid)
  qsubget(jobid{n});
end
%% figure
figure(2)
set(gcf,'position',[50 50 1600 800])
clf
h(3)=plot(snr,ones(size(snr))*(1/100-.5)*-2,'--','color',[1 1 1]*.4,'linewidth',2);
hold on
h(1)=plot(snr,nanmean(skw(:,:,1)),'r','linewidth',2);
plot(snr,quantile(skw(:,:,1),.25,1),'r')
plot(snr,quantile(skw(:,:,1),.75,1),'r')
h(2)=plot(snr,nanmean(skw(:,:,2)),'linewidth',2);
plot(snr,quantile(skw(:,:,2),.25,1))
plot(snr,quantile(skw(:,:,2),.75,1))
legend(h,'PA','SWM','True Skewness','location','southeast')

export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs','Fig2'),'-transparent','-pdf');