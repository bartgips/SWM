function [stats]=bg_bootstrap_interpolate(shapeMat, numIt, frac, verbose, numExtr, fignum, smoothflag)
% stats = bg_bootstrap_interpolate(shapeMat, numIt, frac, verbose, numExtr, fignum, smoothflag)
%
% Uses interpolation together with detection of extrema to calculate skweness index based on T_up/T_down ratio.
% Like other bootstrap functions estimates confidence interval of skewness of noisy shapes contained in
% 'shapeMat' using bootstrap resampling.
%
% This method is useful when the individual shapes are very
% noisy, and this noise hinders you in determining the true skewness of
% every shape individually. Averaging over a large enough (frac) sub sample,
% will average out most of the noise.
%
% Note: this method assumes a unimodal distribution of skewness.
% If the input is a result (a single cluster) from bg_SWM and/or bg_SWM_SA
% this is a fair assumption.
%
% %%%%%%%%%
% input:
% %%%%%%%%%
%
% shapeMat: Matrix of size numShapes x numTimepoints. E.g. the resulting
%           matrices from bg_swm_extract.
%
% numIt:    How many of these resampled distributions should be constructed
%
% frac:     Fraction of total number of samples used every for every bootstrap
%           resample. (i.e. frac*numShapes is the amount of shapes in
%           every subsample)
%           (default = 1; for correct bootstrap statistics this *should be*
%           1)
%
% verbose:  Flag that determines whether the progress (i.e. iteration
%           counter) is written to the terminal/ command window.
%           (default: verbose = 1)
%
% numExtr:  Number of extrema detected to estimate the skewness.
%           (default = 3; should be 3 or greater)
% 
% fignum:   (optional) Figure number/handle to which to plot every sample;
%           smoothed together with the 3 detected extrema.
%           If left empty or set to zero, no figure is plotted.
%
% %%%%%%%%%
% output:
% %%%%%%%%%
%
% stats:  a structure containing the calculated statistics:
% firstly it contains two fields for visualization:
% .meanShape: the mean shape used for detecting asymmetry (interpolated)
% .extrema:   the median of the positions of the extrema used for the calculation
%
% stats contains 2 sub structures .skw and .period. These contain the
% statistics on both skewness and period respectively.
%
% These substructures contain the following fields:"
%
% .mu:    Estimated mean
% .sem:   Estimated standard error of .mu
% .distr: Samples generated by the bootstrap procedure that are used to
%         calculate .mu and .sem.
% .p_t:   p-value for student-t test for rejecting the hypothesis that the
%         mean is equal to zero. (alpha =.05)
% .CI:    Similar to .p_t, but now the 95% confidence interval. If this
%         contains zero, H0 cannot be rejected

if nargin<7
  smoothflag=0;
end

if nargin<3 || isempty(frac)
  frac=1;
end

if frac~= 1
  warning('The sample size for the bootstrap statistics is not equal to the original sample size. Statistics on the mean will not be correct.')
end

if nargin<4
  verbose=1;
end
reverseStr=[];

if nargin >5 && fignum
  figure(fignum)
  clf
  %   set(fignum,'visible','off')
end

if nargin <5
  numExtr=3;
end



[numTemp, tempLen]=size(shapeMat);

skwIdx=nan(numIt,1);
brd=nan(numIt,numExtr);
sampsz=round(frac*numTemp);

% find bias for bg_skewness_pktg_smooth; i.e. find period with least
% variance
nfft=2^nextpow2(tempLen)*4;
ftshape=fft(nanmean(shapeMat),nfft);
[~,midx]=max(abs(ftshape(1:nfft/2+1)));
shapeLen=round(nfft/midx*numExtr/2);

if shapeLen < .65 * size(shapeMat,2)
  varShape=nanvar(shapeMat);
  % push shapes towards centre (increase cost of edges by 25%)
  parabola=[1:numel(varShape)]-numel(varShape)/2-.5;
  parabola=parabola/parabola(end);
  parabola=1+parabola.^2*.25;
  
  bias=conv(varShape.*parabola,ones(1,shapeLen),'valid');
  
  [~,bias]=min(bias);
  bias=[bias bias+shapeLen/2+.5 bias+shapeLen];
else % if number of extrema only just fits in the window, try to centre them
  bias=size(shapeMat,2)/2;
end

interpFac=1e2;
amplitudes=nan(numIt,1);
sigma=0;
for iter=1:numIt
  
  sel=ceil(rand(sampsz,1)*numTemp);
  meanShape=nanmean(shapeMat(sel,:))';
  
  if verbose
    msg=sprintf(['Iteration %d/%d\n'], [iter numIt]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  %% calculating SkwIdx
    if smoothflag
      [skwIdx(iter), brd(iter,:), meanShapeInt, sigma]=bg_skewness_pktg_smooth(meanShape,bias,interpFac,numExtr,sigma);
    else
      [skwIdx(iter), brd(iter,:), meanShapeInt]=bg_skewness_pktg_harsh(meanShape,bias,interpFac,numExtr);
    end
    
    if iter==1
      % make sure that bootstrapping will always focus on the same period
      bias= brd(iter,:)/interpFac;
    end
    
    brddum=brd(iter,:);
    
    amplitudes(iter)=mean(abs(diff(meanShapeInt(brddum(~isnan(brddum))))))/2;
    
    if nargin >5 && fignum
      tInt=[1:numel(meanShapeInt)]/interpFac;
      current_figure(fignum)
      plot(tInt,meanShapeInt)
      vline(brd(iter,:)/interpFac)
      title(sprintf(['Iteration %d/%d\n'], [iter numIt]))
      title(['Iteration ' num2str(iter) '/' num2str(numIt) '; Skewness: ' num2str(skwIdx(iter),'%1.3f')])
      xlim([1 numel(meanShapeInt)]/interpFac)
      drawnow
    end
  
  
  
end



%% skewness
stats.skw.mu=nanmean(skwIdx);
stats.skw.sem=sqrt(numIt/(numIt-1)*nanvar(skwIdx,1));
stats.skw.distr=skwIdx;


% perform t-test on difference from zero (only works if frac=1).
t=stats.skw.mu/stats.skw.sem;
stats.skw.p_t=1-tcdf(abs(t),numTemp-1);
% 95% confidence
alpha=.05;
try
  stats.skw.CI= quantile(skwIdx,[alpha/2 1-alpha/2]);
catch
  warning('no confidence interval(s) calculated.')
end

%% period
periods=diff(brd(:,[1 3]),1,2);
stats.period.mu=nanmean(periods)/interpFac;
stats.period.sem=sqrt(numIt/(numIt-1)*nanvar(periods,1))/interpFac;
stats.period.distr=periods/interpFac;


% perform t-test on difference from zero (only works if frac=1).
t=stats.period.mu/stats.period.sem;
stats.period.p_t=1-tcdf(abs(t),numTemp-1);
% 95% confidence
alpha=.05;
try
  stats.period.CI= quantile(periods,[alpha/2 1-alpha/2]);
end

%% amplitude
stats.amplitude.mu=nanmean(amplitudes);
stats.amplitude.sem=sqrt(numIt/(numIt-1)*nanvar(amplitudes,1));
stats.amplitude.distr=amplitudes;

%% mean extrema positions and shape
stats.extrema=[quantile(brd,.5) ];
meanShapeDum=nanmean(shapeMat)';
t=1:tempLen;
tint=linspace(1,tempLen,numel(meanShapeInt));
stats.meanShape=spline(t',meanShapeDum,tint');
stats.sampSz=sampsz;
% cutout_dum=max(cutout_dum,1);
% cutout_dum=min(cutout_dum,numel(tint));
% stats.meanShape=stats.meanShape(cutout_dum(1):cutout_dum(2));
% stats.extrema=stats.extrema-cutout_dum(1)+1;

function current_figure(h)
set(0,'CurrentFigure',h)


function Y=quantile(x, p)

Y=nan(numel(p),size(x,2));
for k=1:size(x,2);
  xdum=x(:,k);
  xdum=xdum(~isnan(xdum));
  xdum=sort(xdum,1);
  L=size(xdum,1);
  for n=1:numel(p)
    idx=p(n)*(L-.5)+.5;
    remainder=rem(idx,1);
    Y(n,k)=(1-remainder)*xdum(floor(idx))+(remainder)*xdum(floor(idx)+1);
  end
end




