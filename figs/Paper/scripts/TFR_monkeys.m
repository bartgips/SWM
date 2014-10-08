clear all

addpath('/home/mrphys/bargip/OrientationDataMonkeys/matlab_preproc')


basedir={'/home/mrphys/bargip/OrientationDataMonkeys/Okkie';...
          '/home/mrphys/bargip/OrientationDataMonkeys/Spock'};
        
penSel={'009','142'};



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


fdir='/home/mrphys/bargip/GIT/SWM/figs/Paper/tfr';
if ~exist(fdir,'file')
  mkdir(fdir)
end

%%
doTFR=0;
foi=[15:1:80];

for dirit=4;%1:numel(pendirs)
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
    
    chanSelString=cell(1,numel(chansel{dirit}));
    for n=1:numel(chansel{dirit});
      chanSelString{n}=num2str(chansel{dirit}(n),'AD%02d');
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
      cfg.toilim = [-.5 2]; % select 200 ms after stimulus onset;
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
    
    if doTFR
      nandum=isnan(dat);
      dat(nandum)=0;
      clear nandum
      
      % calculate TFRs in chunks of trials such that I don't need the
      % entire torque cluster worth of RAM
      Chunks=1:150:size(dat,1);
      chnLen=ones(size(Chunks))*100;
      chnLen(end)=size(dat,1)-Chunks(end)+1;
      numChunks=numel(Chunks);
      toi=-.3:5/Fs:1.5;
      %%%%% THINK ABOUT KEEPING TRACK OF NaNs %%%%%%%%%%%%%%%
      TFR=zeros(numel(foi),numel(toi));
      for chn=1:numChunks
        [tfr,f,t]=ft_specest_tfr(dat([0:chnLen(chn)-1]+Chunks(chn),:),cfg.toilim(1):1/Fs:cfg.toilim(2),'freqoi',foi,'timeoi',toi,'width',12);
        TFR=TFR+squeeze(nansum(abs(tfr).^2));
      end
      TFR=TFR/size(dat,1);
      
      tfr=bsxfun(@times,TFR,f(:).^2);
      
      figure(99)
      clf
      set(99,'position',[100 100 800 600])
      imagesc(t,f,tfr,[0 nanmax(nanmax(tfr(f>25,t>.2)))]); axis xy;
      xlim([-.3 1])
      vline(.2,'-w','start data for SWM')
      shanknum=ceil(elec/16);
      title(sprintf([pendir{1} '\nContact ' num2str(elec - (shanknum-1)*16) '; shank ' num2str(shanknum)]))
      h=colorbar;
      set(get(h,'ylabel'),'string','power')
      
      xlabel('time rel to stim onset (s)')
      ylabel('frequency')
      
      
      
      fname=[Monkey '_' pendir{1}];
      if elec==1
        export_fig(fullfile(fdir,fname),'-pdf')
      else
        export_fig(fullfile(fdir,fname),'-pdf','-append')
      end
    end
    
    
    
    
  end
  
  
  
end