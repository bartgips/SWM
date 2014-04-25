cd /home/mrphys/bargip/GIT/SWM/figs/stats

%%
skew=.8;
tLen=100;
verbose=0;


snr=[.01:.01:.1 .2:.1:1];
stat=cell(numel(snr),1);

numIt=50;

jobid=cell(numIt,1);
for iter=1:numIt
  jobid{iter}=qsubfeval('bg_sawtooth_shapefind_snr_bat',snr,skew,tLen,1,'memreq',1024^2*100,'timreq',60*200);
% [skw(iter,n,:)]=bg_sawtooth_shapefind_snr(snr(n),skew,tLen,verbose);

end

save jobid_8skew jobid snr

%%
load jobid

%%%%%change this

skw=nan([numel(snr) numel(jobid) 2]);
%%%%%%%
missed=[];
for n=1:numel(jobid)
  home
  disp(['loading result ' num2str([n numel(jobid)],'%d/%d')])
  try
    argout=bg_qsubget(jobid{n});
    skw(:,n,:)=argout{1};
  catch
    missed=[missed n];
  end
end

save('run10s_.01_1.mat','skw','missed','snr')

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
export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs','Fig2'),'-transparent','-pdf');