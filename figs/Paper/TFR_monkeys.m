clear all

addpath('/home/mrphys/bargip/OrientationDataMonkeys/matlab_preproc')


basedir={'/home/mrphys/bargip/OrientationDataMonkeys/Okkie';...
  '/home/mrphys/bargip/OrientationDataMonkeys/Spock'};

penSel={'009','142'};

Monkey={'Monkey O','Monkey S'};

chansel={1:48, 17:32};

alignment=  [9 2 3 3; 141 0 2 2];
%% okkie Chansel

if strcmp(basedir(end-4:end),'Okkie')
  chansel={1:16
    1:48
    1:48
    1:48
    1:16};
  %1:16};
  Fhp=25;
  
  Monkey='Okkie';
  savdirbase='/home/mrphys/bargip/OrientationDataMonkeys/Torque/Okkie';
  alignment=  [7 3 2 3
    9 2 3 3
    38 2 3 5
    53 1 2 4];
elseif strcmp(basedir(end-4:end),'Spock')
  chansel={17:32
    17:32
    17:32
    1:16
    1:16
    1:16
    1:48
    1:48
    1:48
    1:48};
  Fhp=24;
  Monkey='Spock';
  
  savdirbase='/home/mrphys/bargip/OrientationDataMonkeys/Torque/Spock';
  alignment=  [160 5 4 4
    161 5 4 3
    162 3 3 3
    139 0 NaN 8
    141 0 2 2
    142 0 NaN 1
    143 1 7 NaN
    144 2 5 5];
end
%%

fdir='/home/mrphys/bargip/GIT/SWM/figs/Paper/tfr';
if ~exist(fdir,'file')
  mkdir(fdir)
end

%%
foi=[5:1:65];
TFR=cell(numel(basedir),1);

%%
for monkeyIt=1:numel(basedir)
  pendir=dir(basedir{monkeyIt});
  pendir={pendir.name};
  pendir=bg_selectCell(pendir,'pen','~.mat');
  pendir=bg_selectCell(pendir,penSel{monkeyIt});
  fnames=dir(fullfile(basedir{monkeyIt},pendir{1}));
  fnames={fnames.name};
  plexon_data=bg_selectCell(fnames,'.nex');
  
  recordn=bg_selectCell(fnames,'record');
  
  
  
  data_lfp=cell(24,numel(plexon_data));
  
  for k=1:numel(plexon_data)
    recordn2=fullfile(basedir{monkeyIt},pendir{1},recordn{1});
    plexon_data2=fullfile(basedir{monkeyIt},pendir{1},plexon_data{1});
    fnames=dir(recordn2);
    expdata=bg_selectCell({fnames.name},'.mat');
    expdata=fullfile(basedir{monkeyIt},pendir{1},recordn{1},expdata{1});
    
    load(expdata)
    
    %   output_dir=cd;
    trl_data=exp_data.trial_data{1};
    ctx_data=trl_data{1};
    ctx_codes=ctx_data(:,1);
    ctx_times=ctx_data(:,2);
    sample_ctx_codes=ctx_codes(1:4);
    
    exp_pars=exp_data.exp_info;
    if ~isfield(exp_pars,'post_stim_time')
      exp_pars.post_stim_time=0;
    end
    
    ctx_triggers=exp_pars.triggers;
    clear trl_data
    
    %%%%%%%%% read in the plexon_data and go though the trials
    
    cfg.dataset=(plexon_data2);
    
    cfg.sample_ctx_codes=sample_ctx_codes;
    %%%%% read the header%%%%%%%%%%%%%%%%
    
    hdr = ft_read_header(cfg.dataset);
    
    
    
    ch_labels=hdr.label; %find a way to select the spike channels from this
    
    
    
    LFP_channels=cat(1,ch_labels{~cellfun(@isempty,strfind(ch_labels,'AD'))});
    %   sp_channel_index=find(~cellfun(@isempty,regexp(ch_labels,'sig0\S+_wf')));
    %   sp_channels=cat(1,ch_labels{sp_channel_index});
    %   spk_channel_index=[sp_channel_index str2num(sp_channels(:,4:6))];
    
    
    if ~isfield(exp_pars,'V2_electrodes')
      exp_pars.V2_electrodes=[];
    end
    
    electrode_names=str2num(LFP_channels(:,3:end));
    
    
    % script to find the number of conditions
    [~, col]= find(exp_data.ndir);
    find_cnds= col;
    
    clear exp_data
    
    chanSelString=cell(1,numel(chansel{monkeyIt}));
    for n=1:numel(chansel{monkeyIt});
      chanSelString{n}=num2str(chansel{monkeyIt}(n),'AD%02d');
    end
    
    %%%%% get the trials for this condition
    for n=1:24
      cfg = [];
      cfg.dataset=[plexon_data2];
      cfg.condition = n;%find_cnds([8,16]); % conidtions 8 and 16 ( 45 and 135)
      cfg.sample_ctx_codes=sample_ctx_codes;
      cfg.trialfun = 'get_trigger_data';
      cfg.ctx_triggers=ctx_triggers;
      cfg.exp_pars=exp_pars;
      cfg.exp_pars.baseline= -2;  % how much baseline
      cfg.whole_trial=1;
      cfg = ft_definetrial(cfg);
      
      
      
      cfg.hpfilter='no';
      %       cfg.hpfreq=Fhp; % HP filter to remove alpha
      cfg.channel= chanSelString;
      %   cfg.dftfilter = 'yes'; % optional
      if ~isempty(cfg.trl)% && ((V1_channels(1,slot)==0) || (V2_channels(1,slot)==0))
        data_lfpdum = ft_preprocessing(cfg);
      end
      
      cfg=[];
      cfg.toilim = [-1 2]; % select 200 ms after stimulus onset;
      data_lfp{n,k} = ft_redefinetrial(cfg, data_lfpdum);
      clear data_lfpdum
    end
  end
  
  % concatenating all the data
  % determining the size of the matrix (numtimpoints is actually bounded to
  % 1801 right now, because of cfg.toilim above. So an unnecessary step, but this is more flexible)
  numtrls=0;
  numtimpoints=0;
  for q=1:numel(data_lfp)
    numtrls=numtrls+numel(data_lfp{q}.trial);
    for qq=1:numel(data_lfp{q}.trial)
      numtimpoints=max([numtimpoints size(data_lfp{q}.trial{qq},2)]);
    end
  end
  
  clear ch_labels ctx_triggers exp_pars
  Fs=hdr.Fs;
  clear hdr
  
  
  numelec=numel(chanSelString);
  numShanks=ceil(numelec/16);
  penNo=str2double(pendir{1}(end-2:end));
  addpath /home/mrphys/bargip/test/Matlab/
  if any(alignment(:,1)==penNo)
    aligOff=alignment(alignment(:,1)==penNo,2:end);
  else
    aligOff=0;
  end
  
  
  
  dat=nan(numtrls*16,numtimpoints);
  cnt=0;
  for q=1:numel(data_lfp)
    for qq=1:numel(data_lfp{q}.trial)
      for qqq=1:16
        cnt=cnt+1;
        dat(cnt,1:numel(data_lfp{q}.trial{qq}(qqq,:)))= data_lfp{q}.trial{qq}(qqq,:);
      end
    end
  end
  
%   dat=bg_filter_brickwall(dat.', Fs, [51 49]).';
  
    nandum=isnan(dat);
    dat(nandum)=0;
    
    % calculate TFRs in chunks of trials such that I don't need the
    % entire torque cluster worth of RAM
    chunkSize=150;
    chunks=1:chunkSize:size(dat,1);
    chnLen=ones(size(chunks))*chunkSize;
    chnLen(end)=size(dat,1)-chunks(end)+1;
    numChunks=numel(chunks);
    toi=-.3:5/Fs:1.5;
    tAxis=[cfg.toilim(1):1/Fs:cfg.toilim(2)];
    TFR{monkeyIt}=zeros(numel(foi),numel(tAxis));
    dat=dat(randperm(size(dat,1)),:); % useful when only using a subset of the data
    for chn=1:10%numChunks
      [tfr,f,t]=ft_specest_tfr(dat([0:chnLen(chn)-1]+chunks(chn),:)/1e3,tAxis,'freqoi',foi,'timeoi',tAxis,'width',12);
      %  controlling for missing data during averaging
      nandum2=sum(~nandum([0:chnLen(chn)-1]+chunks(chn),:))/chnLen(chn);
      addDum=squeeze(nansum(abs(tfr).^2));
%       addDum=bsxfun(@times,addDum,1./nandum2);
%       addDum(isnan(addDum))=0;
      TFR{monkeyIt}=TFR{monkeyIt}+addDum;
    end
    nandum3=all(isnan(tfr),1);
    TFR{monkeyIt}=TFR{monkeyIt}/min(sum(chnLen(1:chn)),size(dat,1));
    TFR{monkeyIt}(nandum3)=nan;
    
    
    
    

  
  
  
  
end
%%
save TFRs TFR t f
%%

load /home/mrphys/bargip/GIT/SWM/figs/Paper/TFRs.mat

%%
figure(99)
clf
set(99,'position',[100 100 1400 600])
fsz=14;
sh=.09;
for n=1:numel(TFR)
  subaxis(1,2,n,'sh',sh)
  plotDum=TFR{n};
  plotDum=bsxfun(@times,plotDum,f(:).^2);
%   plotDum=(bsxfun(@rdivide,plotDum,nanmean(plotDum(:,t<-.1),2)));
  plotDum=(bsxfun(@minus,plotDum,nanmean(plotDum(:,t<-.1),2)));
  imagesc(t,f,(plotDum),[0 max(max(plotDum(f<55,:)))])
  set(gca,'fontsize',fsz)
  [hl, ht]=vline(-.1,'--w','baseline');
  set(hl,'linewidth',2,'color',[1 1 1]*.99);
  txtpos=get(ht,'position');
  txtpos(1)=-.4;
  set(ht,'position',txtpos,'fontsize',fsz)
  axis xy;
  xlim([-.5 1])
%   vline(.2,'--w','start data for SWM')
  title(Monkey{n},'fontsize',fsz)
  bg_figureLabel(char(64+n))
  h=colorbar;
  set(h,'fontsize',fsz)
  set(get(h,'ylabel'),'string','LFP power (V^2\cdots^2)','fontsize',fsz)
  xlabel('time relative to stimulus onset (s)','fontsize',fsz)
  ylabel('frequency (Hz)','fontsize',fsz)
end
colormap hot


%%
export_fig(fullfile('/home/mrphys/bargip/GIT/SWM/figs/Paper','Fig4_monkeyTFRs'),'-transparent','-png','-eps','-pdf',99)