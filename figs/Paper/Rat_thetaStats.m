%% compiling/loading function
compileBS=0;
if compileBS
  temp=pwd;
  cd ~/GIT/SWM/related_functions/compiled/
  revNum=importdata('/home/mrphys/bargip/GIT/SWM/version');
  revNum=regexp(revNum,'(r\S+-\d+)-','tokens');
  revNum=strrep(revNum{1}{1}{1},'-','_');
  batchid=['bg_bootstrap_interpolate_' revNum];
  compiledfun=qsubcompile({@bg_bootstrap_interpolate},'toolbox',{'stats'},'batchid',batchid);
  compiledfun.revNum=revNum;
  fname=fullfile('~/GIT/SWM/related_functions/compiled',batchid);
  save(fname,'compiledfun');
  save('~/GIT/SWM/related_functions/compiled/latest_bootstrap','compiledfun')
  cd(temp)
else
  load ~/GIT/SWM/related_functions/compiled/latest_bootstrap
  revNum=compiledfun.revNum;
end

%% loading SWM results
basedir='/home/mrphys/bargip/Data_Eric/Raw_Data/Colgin_Hippocampus';
cfgs=cell(3,2);
for n=1:3
  for k=1:2
    fname=fullfile(basedir,'Batch','2014_08','output',['cfg_' num2str([n k],'%d_%d')]);
    load(fname)
    cfgs{n,k}=cfg;
  end
end

%%
label=cell(1,numel(cfgs));
% meanshape=nan(cfgs{1}.winLen,numel(cfgs));
stats=cell(size(cfgs));
verbose=0;
fignum=[];
numIt=1000;
frac=1;
batchid=['RatsAsymmetryBS_' revNum];
jobid=cell(size(cfgs));

for k=1:numel(cfgs)
  
  cfg=cfgs{k};
  
  
  labeldum=cfg.fname;
  [labeldum,fnamedum]=fileparts(labeldum);
  [~,labeldum]=fileparts(labeldum);
  label{k}=labeldum(1:end-3);
  
  cfg.newLen=round(1.2*cfg.winLen);
    
  [z, s]=bg_swm_extract(cfg);
%   z=z{1};
  zsel=~any(isnan(z),2);
  cost=sum(bsxfun(@minus,z,nanmean(z)).^2,2);
  zsel=cost<median(cost);
  
  
  
  stats{k}=bg_bootstrap_interpolate(z(zsel,:),numIt, frac, 1, 4);
%   jobid{k}=qsubfeval(compiledfun, z(zsel,:), numIt, frac, verbose, 'memreq', 1.5*1024^3, 'timreq', 10*3600, 'batchid', batchid);
%   stats{k}.meanShape=nanmean(z(zsel,:));
  stats{k}.label=label{k};

end

cd /home/mrphys/bargip/Data_Eric/Raw_Data/Colgin_Hippocampus/Batch/2014_08/stats
save(['statsDum_' revNum ],'stats','jobid')

%%
stats=bg_concat_struct(stats);

%% plotting

col=lines(3);
xl=.6;
lw=2;

fs=2e3;
legLabels={'Rat1', 'Rat2', 'Rat3'};
tVec=[1:size([stats.meanShape],1)]/100/fs*1e3-20;

figure(27)
clf
set(gcf,'position',[100 100 600 800])

subselstat=1:3;
skw=bg_concat_struct(stats.skw);
skw=reshape(skw,[3,2]);
dum=[skw.distr];
bins=minmax(dum(dum<1));
% bins=[-.5 0];
bins=linspace(bins(1),bins(2),1e2);

subaxis(3,1,1,'ML',.08,'MR',.08,'sv',.08)
hold on
h=nan(1,2);
for n=subselstat
h(n)=plot(tVec,[stats(n,1).meanShape],'linewidth',lw,'color',col(n,:));
plot(tVec(round(stats(n,1).extrema)),stats(n,1).meanShape(round(stats(n,1).extrema)),'o','color',col(n,:),'markersize',10)
end
title('Hippocampal Gyrus','fontsize',14)
ylabel('averaged z-scored LFP')
xlabel('time (ms)')
axis tight
ax=axis;
yrng=ax(3:4);
ylim(yrng+[-.1 .1]*diff(yrng));
legend(h,{stats(:,1).label},'location','southwest')
legend(h,legLabels,'location','best')
hh=vline([0 200],'--k');
set(hh,'linewidth',2)
bg_figureLabel('A')

subaxis(3,1,2)
hold on
for n=subselstat
h(n)=plot(tVec,[stats(n,2).meanShape],'linewidth',lw,'color',col(n,:));
plot(tVec(round(stats(n,2).extrema)),stats(n,2).meanShape(round(stats(n,2).extrema)),'o','color',col(n,:),'markersize',10)
plot(tVec(round(stats(n,2).extrema)),stats(n,2).meanShape(round(stats(n,2).extrema)),'x','color',col(n,:),'markersize',10)
end

ylabel('averaged z-scored LFP')
title('Area CA1','fontsize',14)
xlabel('time (ms)')
axis tight
ax=axis;
yrng=ax(3:4);
ylim(yrng+[-.1 .1]*diff(yrng));
% legend(h,{stats(:,1).label},'location','southwest')
hh=vline([0 200],'--k');
set(hh,'linewidth',2)
bg_figureLabel('B')

subaxis(3,2,1,3)
h=barplotsemstar([skw(subselstat,1).mu],[skw(subselstat,1).sem],[skw(subselstat,1).p_t]<.05,col);
ylabel('Skewness Index')
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',{stats(:,1).label})
set(gca,'xticklabel',legLabels)
title('Hippocampal Gyrus','fontsize',12)
bg_figureLabel('C')
% subaxis(2,4,2,2)
% N=nan(numel(bins),2);
% for n=subselstat
%   N(:,n)=hist(skw(n).distr,bins);
% end
% plot(bins,N/sum(N(:)),'linewidth',2)
% xlim(minmax(bins))
% xlabel('skewness')




subaxis(3,2,2,3)
h=barplotsemstar([skw(subselstat,2).mu],[skw(subselstat,2).sem],[skw(subselstat,2).p_t]<.05,col);
% ylabel('Skewness Index')
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',{stats(:,2).label})
set(gca,'xticklabel',legLabels)
title('Area CA1','fontsize',12)
% subaxis(2,4,4,2)
% N=nan(numel(bins),2);
% for n=subselstat
%   N(:,n)=hist(skw(subselstat(n),2).distr,bins);
% end
% plot(bins,N/sum(N(:)),'linewidth',2)
% xlim(minmax(bins))
% xlabel('skewness')

% mtit(strrep(['skew_Bootstrap_' num2str(frac) 'f_' num2str(numel(stats(1).skw.distr)) ],'_','\_'),'yoff',.05)
bg_figureLabel('D')
savefigs=0;
if savefigs
  fname=['Fig3_RatStats_' revNum];
  fdir='/home/mrphys/bargip/GIT/SWM/figs/Paper';
  export_fig(fullfile(fdir,fname),'-transparent','-png','-pdf','-eps')
  
end

%% TFR
[TFR, sd, sem, PLV]=deal(cell(size(cfgs)));
for n=1:numel(cfgs)
cfg=cfgs{n};
cfg.newLen=cfg.winLen*2;





[TFR{n}, f, sd{n}, sem{n}, PLV{n}]=bg_swm_TFR(cfg, 20:5:200);


end

%%
figure (1); 
clf
set(gcf,'position',[100 100 1200 600])
labels={'Rat 1', 'Rat 2', 'Rat 3'; 'Rat 1', 'Rat 2', 'Rat 3'}';
for n=1:numel(TFR)
  [a,b]=ind2sub(size(TFR),n);
  TFRnorm=bsxfun(@minus,TFR{n},mean(TFR{n},2));
TFR_t=TFRnorm./sem{n};
TFR_p=1-tcdf(abs(TFR_t),cfg.numTemplates-1);
TFR_t_mask=TFR_t;
mask=(TFR_p<.05/(numel(TFR_t)));
  cfg=cfgs{n};
cfg.newLen=cfg.winLen*2;

  [z,s,cfg]=bg_swm_extract(cfg);
  t=[1:cfg.newLen]-cfg.addLen;
t=t/cfg.fs;
subaxis(2,3,n,'sv',.08)
imagesc(t,f,TFR_t_mask,maxabs(TFR_t)); 
set(gca,'ydir','normal')
hold on
plot(t,nanmean(z)*std(f)+mean(f),'k','linewidth',2)
% [~,tt]=fileparts(cfg.fname);
% title(tt)
if b==1
title(labels{n},'fontsize',2)
end
h=colorbar;
set(get(h,'ylabel'),'string','t-value');
aa=(get(gca,'children'));
set(aa(2),'alphaData',mask)
h=vline([0 cfg.winLen/cfg.fs],'--k');
set(h,'linewidth',2);
if b==2
xlabel( 'time (ms) ')
end
end

axes('position',[0,0,1,1],'visible','off');
h=text(.05,.6,'Hippocampal Gyrus','rotation',90,'fontsize',16);
h=text(.05,.25,'Area CA1','rotation',90,'fontsize',16);

savefigs=0;
if savefigs
  fname=['RatTFR_' revNum];
  fdir='/home/mrphys/bargip/Data_Eric/Raw_Data/Colgin_Hippocampus/Batch/2014_08/figs';
  export_fig(fullfile(fdir,fname),'-transparent','-png','-pdf','-eps',1)
  
end