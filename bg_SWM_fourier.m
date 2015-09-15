function cfg=bg_SWM_fourier(cfg,dat)
% cfg=bg_SWM_fourier(cfg,dat)
% 
% emulates the functionality of bg_SWM.
% instead of the rigorous SWM algorithm, it uses band-pass filtering to
% detect 0-crossings in the phase evolution.
% i.e. it centres windows on cosine peaks.

if nargin<2
  dum=load(cfg.fName, cfg.varName);
  eval(['dat=dum.' cfg.varName ';'])
end

if iscolumn(dat);
  dat=dat.';
end

if isfield(cfg,'foi')
  foi=cfg.foi;
else
  error('no frequency given for bp filtering (cfg.foi)')
end

if isfield(cfg,'bw')
  bw=cfg.bw;
else
  warning('no bandwidth given for bp filtering (cfg.bw), using +-2Hz')
  bw=2;
end

if isfield(cfg,'fs')
  fs=cfg.fs;
else
  error('no sampling frequency given (cfg.fs)')
end

if isfield(cfg,'winLen')
  winLen=cfg.winLen;
else
  winLen=fs*1/foi;
  cfg.winLen=winLen;
end

if isfield(cfg,'type')
  type=cfg.type;
else
  type='gauss';
  cfg.type=type;
end

if isfield(cfg,'winPerTrial')
  cfg.maxlocs=cfg.winPerTrial;
end

if isfield(cfg,'MEGflag')
  % MEG-signals have "arbritrary" polarity, therefore they can be rectified 
  % by flipping phases 180 degrees if they are in antiphase with the 
  % strongest signal.
  MEGflag=cfg.MEGflag;
else
  MEGflag=false;
end


switch type
  case 'gauss'
    dt = 1/fs;
    sf = bw;
    st = 1/(2*pi*sf);
    t = -5*st:dt:5*st;
    A=1/(sqrt(2*pi)*st)*2/fs;
    filtkern = A*(exp(-t.^2/(2*st^2)))';
    filtkern = filtkern .* exp(-1i*t'*2*pi*foi);
    analytic=convn(dat,filtkern','same');
  case 'butt'
    filtdat=ft_preproc_bandpassfilter(dat,fs,[-bw bw]+foi);
    analytic=hilbert(filtdat.').';
  case 'brick'
    if numel(foi)~=2
      error('cfg.foi should contain 2 frequencies when using brickwall filter')
    end
    % brickwall complex BP filters
    L=size(dat,2);
    df=1/L*fs;
    if mod(L,2)
      f=[0:df:fs/2 -fs/2:df:0-df]';
    else
      f=[0:df:fs/2-df -fs/2:df:0-df]';
    end    
    filtFT=zeros(numel(f),1);
    filtFT(f>=foi(1) & f<=foi(2))=1;
    analytic=ifft(bsxfun(@times,fft(dat.'),filtFT)).';
end



if MEGflag && ~ismatrix(analytic)
  %find strongest signal (most power)
  [~,mainSens]=max(squeeze(sum(sum(abs(analytic).^2,2))));
  cohDum=bsxfun(@times,analytic,conj(analytic(:,:,mainSens)));
  cohDum=squeeze(sum(sum(cohDum.*abs(cohDum)),2)); % weight by power
  signFlip=sign(real(cohDum));
  
  analytic(:,:,:)=bsxfun(@times,analytic(:,:,:),permute(signFlip,[3 2 1]));
end

if ~ismatrix(analytic)
  analytic=sum(analytic(:,:,:).*abs(analytic(:,:,:)),3); %weight by power, rather than amplitude
end
 
%finding oscillation peaks;
ph_tot=unwrap(permute(angle(analytic),[2 1 3:ndims(analytic)]));
numlocs=max(floor(max(ph_tot(:,:))/(2*pi)));
numtrials=size(ph_tot,2);
loc=nan(numtrials,numlocs);


for trialsel=1:numtrials;
  ph=squeeze(ph_tot(:,trialsel));
  peaks=(2*pi):(2*pi):max(ph);
  for peakidx=1:numel(peaks)
    [~, loc(trialsel,peakidx)]=min(abs(ph-peaks(peakidx)));
  end
end

if isfield(cfg,'maxlocs') %only extract highest power
  maxlocs=cfg.maxlocs;
  cfg.best_loc=nan(numtrials,maxlocs);
  for trialsel=1:numtrials;
    nansel=~isnan(loc(trialsel,:));
    locdum=loc(trialsel,nansel);
    amp=abs(analytic(trialsel,locdum));
    [a, sel]= sort(amp,2,'descend');
    sel=sel(~isnan(a));
    if numel(sel)>maxlocs
      sel=sel(1:maxlocs);
    end
    cfg.best_loc(trialsel,1:numel(sel))=sort(locdum(sel));
  end
  cfg.best_loc=cfg.best_loc-floor(winLen/2);
  cfg.best_loc=[cfg.best_loc ones(numtrials,1)*size(dat,2)];
  cfg.bestclust{1}.numWindows=numtrials*maxlocs;
  
else
  cfg.best_loc=loc;
  cfg.best_loc=cfg.best_loc-floor(winLen/2);
  cfg.bestclust{1}.numWindows=numel(loc)-size(loc,1);
end

cfg.bestclust{1}.trl=repmat([1:numtrials]',1,size(cfg.best_loc,2)-1);
cfg.bestclust{1}.tidx=repmat(1:size(cfg.best_loc,2)-1,size(cfg.best_loc,1),1);
cfg.numWindows=numtrials*(size(cfg.best_loc,2)-1);
cfg.bestclust{1}.linidx=1:cfg.numWindows;


    
