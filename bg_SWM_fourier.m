function cfg=bg_SWM_fourier(cfg,dat)
% cfg=bg_SWM_fourier(cfg,dat)
% 
% emulates the functionality of bg_SWM.
% instead of the rigorous SMW algorithm, it uses band-pass filtering to
% detect 0-crossings in the phase evolution.
% i.e. it centres windows on cosine peaks.

if nargin<2
  dum=load(cfg.fname, cfg.varname);
  eval(['dat=dum.' cfg.varname ';'])
end

if isfield(cfg,'fbp')
  fbp=cfg.fbp;
else
  error('no frequency given for bp filtering (cfg.fbp)')
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
  winLen=fs*1/fbp;
  cfg.winLen=winLen;
end

if isfield(cfg,'type')
  type=cfg.type;
else
  type='gauss';
  cfg.type=type;
end

if isfield(cfg,'numTemplates')
  cfg.maxlocs=cfg.numTemplates;
end


switch type
  case 'gauss'
    dt = 1/fs;
    sf = bw;
    st = 1/(2*pi*sf);
    t = -5*st:dt:5*st;
    A=1/(sqrt(2*pi)*st)*2/fs;
    filtkern = A*(exp(-t.^2/(2*st^2)))';
    filtkern = filtkern .* exp(-1i*t'*2*pi*fbp);
    analytic=conv2(dat,filtkern','same');
  case 'butt'
    filtdat=ft_preproc_bandpassfilter(dat,fs,[-bw bw]+fbp);
    analytic=hilbert(filtdat.').';
end


 
%finding alpha peaks;
ph_tot=unwrap(angle(analytic).');
numlocs=max(floor(max(ph_tot)/(2*pi)));
numtrials=size(ph_tot,2);
loc=nan(numtrials,numlocs);


for trialsel=1:numtrials;
  ph=ph_tot(:,trialsel);
  peaks=(2*pi):(2*pi):max(ph);
  %   if min(ph)<0
  %     peaks=[0 peaks];
  %   end
  peaktims=nan(1,numel(peaks));
  peaktimidx=nan(1,numel(peaks));
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
  cfg.bestclust{1}.numtemplates=numtrials*maxlocs;
  
else
  cfg.best_loc=loc;
  cfg.best_loc=cfg.best_loc-floor(winLen/2);
  cfg.bestclust{1}.numtemplates=numel(loc)-size(loc,1);
end

cfg.bestclust{1}.trl=repmat([1:numtrials]',1,size(cfg.best_loc,2)-1);
cfg.bestclust{1}.tidx=repmat(1:size(cfg.best_loc,2)-1,size(cfg.best_loc,1),1);
cfg.numtemplates=numtrials*(size(cfg.best_loc,2)-1);
cfg.bestclust{1}.linidx=1:cfg.numtemplates;


    
