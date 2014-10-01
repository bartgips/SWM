function [skwIdx, stats]=bg_sawtooth_shapefind_snr(snr,skew,tLen,verbose,cfg)
% [skwIdx, stats]=bg_sawtooth_shapefind_snr(snr,skew,tLen,verbose,cfg)
%
% Bart Gips; April 2014

if nargin<4
  verbose=false;
end

if numel(verbose)<2
  verbose=[verbose, verbose];
end

% generating synthetic data to compare wavelet convolution to sliding
% window matching

fs=1e3;
dt=1/fs;
if nargin<3
  tLen=100;
end
t=0:dt:(tLen-dt);
frq=10;

numT=numel(t);

if verbose(1)
  disp('generating signals...')
end


% constant sawtooth;
if nargin<2
  skew=1;
else
  %transform to skewness that sawtooth function uses
  skew=(skew+1)/2;
end
signal=sawtooth(t'*2*pi*frq,skew);

if snr<inf
var_sig=var(signal);

f=linspace(0,fs/2,round(numT/2))';
noisefft=1./(f)*1e3/snr;
phase=rand(size(noisefft))*2*pi;

noisefft=noisefft.*exp(1i*phase);

noisefft(abs(f)<1)=0;
noisefftfull=zeros(numT,1);
noisefftfull(1:numel(f))=noisefft;
noisefftfull(numel(f)+1:end)=conj(flipud(noisefft));
noise=ifft(noisefftfull,'symmetric');

var_noise=var(noise);

noise=sqrt(var_sig/(var_noise*snr))*noise;
else
  noise=0;
end
signalOrig=signal+noise;
% highpass
signal=ft_preproc_highpassfilter(signalOrig',fs,frq*.75,5);

%% SWM

if verbose(1)
  disp('Sliding window matching...')
end

winLen=fs/frq*2;
% guard=round(.75*winLen);
[meanShapeZ, meanStdZ]=deal(nan(winLen,2));

if nargin<5
  cfg=[];
end
cfg.winLen=winLen;
cfg.fs=fs;
% cfg.guard=guard;
cfg.numIt=2e4;
cfg.fullOutput=1;
cfg.Tfac=logspace(-4,0,20);
cfg.konstant=1e3;
cfg.verbose=verbose(2);
% cfg.Fhp=5;
if numel(verbose)>2
  cfg.dispPlot=verbose(3);
else
  cfg.dispPlot=0;
end
% [cfgSWM]=bg_SWM_SA(cfg,signal);
[cfgSWM]=bg_SWM(cfg,signal);

[z,s]=bg_swm_extract(cfgSWM,signal);

% remove ill-fitters
costdistr=sum(cfgSWM.costDistr{1},2);
remsel=[costdistr>(mean(costdistr)+2*std(costdistr))]*0;

meanShapeZ(:,2)=nanmean(z(~remsel,:));
meanStdZ(:,2)=nanstd(z(~remsel,:));

meanShapeS(:,2)=nanmean(s(~remsel,:));
meanStdS(:,2)=nanstd(s(~remsel,:));


%% fourier method

if verbose(1)
  disp('Alpha phase aligned averaging...')
end

cfg=[];
cfg.fbp=frq;
cfg.bw=2;
cfg.fs=fs;
cfg.winLen=winLen;
cfg.maxlocs=cfgSWM.numTemplates;
cfg.verbose=verbose;
cfg.numTemplates=sum(~remsel);
cfg2=bg_SWM_fourier(cfg,signal);

[z,s]=bg_swm_extract(cfg2,signal);

meanShapeZ(:,1)=nanmean(z);
meanStdZ(:,1)=nanstd(z);

meanShapeS(:,1)=nanmean(s);
meanStdS(:,1)=nanstd(s);

%% calculating skewness
if verbose(1)
  disp('Calculating skewness...')
end

[skwIdx, brd, xout]=bg_skewness_pktg_harsh(meanShapeZ);


stats=[];
% stats.SkwIdx=SkwIdx;
stats.meanShapeZ=meanShapeZ;
stats.meanStdZ=meanStdZ;
stats.meanShape=meanShapeS;
stats.meanStd=meanStdS;
% stats.powspec=pow;
% stats.f=ffrq;
stats.signal=signal;
stats.signalOrig=signalOrig;
stats.t=t;
stats.cfgSWM=cfgSWM;
stats.snr=snr;
stats.skwIdx=skwIdx;
stats.pktg=brd;
stats.meanShapeInt=xout;

