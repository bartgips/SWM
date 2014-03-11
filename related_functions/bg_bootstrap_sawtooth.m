function [stats]=bg_bootstrap_sawtooth(shapeMat, frac, numIt, verbose, fignum)
% stats = bg_jackknife_sawtooth(shapeMat, frac, numIt, verbose, fignum)
%
% Estimates confidence interval of skewness of noisy shapes contained in
% 'shapeMat' using bootstrap resampling.
%
% This method is useful when the individual shapes are very
% noisy, and this noise hinders you in determining the true skewness of
% every shape individually. Averaging over a larg enough (frac) sub sample,
% will average out most of the noise.
%
% Note: this method assumes a unimodal distribution of skewness.
% If the input is a result (a single cluster) from bg_SWM this is a
% fair assumption.
%
% %%%%%%%%%
% input:
% %%%%%%%%%
%
% shapeMat: Matrix of size numShapes x numTimepoints. E.g. the resulting
%           matrices from bg_swm_extract.
% 
% frac:     Fraction of total number of samples used every for every bootstrap
%           resample. (i.e. frac*numShapes is the amount of shapes in
%           every subsample)
% 
% numIt:    How many of these resampled distributions should be constructed
%
% verbose:  Flag that determines whether the progress (i.e. iteration
%           counter) is written to the terminal/ command window. 
%           (default: verbose = 1)
% 
% 
% fignum:   (optional) Figure number/handle to which to plot every sample with
%           its respective sawtooth fit. If left empty or set to zero, no
%           figure is plotted.
%
% %%%%%%%%%
% output:
% %%%%%%%%%
%
% stats:  a structure containing the calculated statistics:
%
% .mu:    Mean skewness
% .sem:   Estimated standard error of .mu
% .distr: Samples generated by the jackknife procedure that are used to
%         calculate .mu and .sem.

if nargin<4
  verbose=1;
end
reverseStr=[];

if nargin >4 && fignum
  figure(fignum)
  clf
  set(fignum,'visible','off')
end

[numTemp, tempLen]=size(shapeMat);
skwIdx=nan(numIt,1);
sampsz=round(frac*numTemp);

%% cut off sides of found shape (focus on centre of window +- 5%)
% (for better fitting
meanShapedum=nanmean(shapeMat)';
meanShapedum=meanShapedum-nanmean(meanShapedum);
% detect dominant frequeny
spec=fft(padarray(meanShapedum,5*tempLen,0,'both'));
spec=spec(1:ceil(end/2));
[~,midx]=max(abs(spec));
shapeLen=11*tempLen/midx;

% shift of shape centre
shapeCentre=round(tempLen/2-atan(imag(spec(midx))/real(spec(midx)))/(2*pi)*shapeLen);

brd=shapeCentre+[-1 1]*shapeLen*.55;
brd=[max(floor(brd(1)),1) min(ceil(brd(2)),tempLen)];

%resize and detrend
shapeMat=shapeMat(:,brd(1):brd(2));
shapeMat=detrend(shapeMat')';

%%
% intial fit parameters
X0=[diff(minmax(nanmean(shapeMat)))/2, tempLen/2, 0, mean(shapeMat(:)), 0]';
fitOptions=optimset('Algorithm','interior-point','Display','off');

for iter=1:numIt
  
  sel=ceil(rand(sampsz,1)*numTemp);
  meanShape=nanmean(shapeMat(sel,:))';
  
  if verbose
  msg=sprintf(['Iteration %d/%d\n'], [iter numIt]);
  fprintf([reverseStr, msg]);
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  %% calculating SkwIdx
  
  %sawtoothfit
  lowerBound=[0; (tempLen/4); -(tempLen/2); -inf; -.999];
  upperBound=[10; (tempLen*3/4); (tempLen/2); inf; .999];
  [X0]=fmincon(@(x)sawtoothfit(meanShape,x),X0,[],[],[],[],lowerBound,upperBound,[],fitOptions);
  skwIdx(iter)=X0(5);
  
  if nargin >4 && fignum
    stdum=X0(1)*sawtooth(([1:numel(meanShape)]-X0(3))*2*pi/X0(2),(X0(5)+1)/2)+X0(4);
    figure(fignum)
    plot(meanShape)
    hold on
    plot(stdum,'r')
    hold off
    title(sprintf(['Iteration %d/%d\n'], [iter numIt]))
    xlim([1 numel(meanShape)])
  end
  
end

stats.mu=mean(skwIdx);
stats.sem=sqrt(numIt/(numIt-1)*var(skwIdx,1));
stats.distr=skwIdx;

function c=sawtoothfit(dat, lambda)
% c=sawtoothfit(dat, lambda)
% 
% calculated cost function of sawtooth with parameters 'lambda' compared to
% data.
% 
% lambda = [ amplitude; period; offset_x; offset_y; skewness];

dat=dat(:);
swth=lambda(1)*sawtooth(([1:numel(dat)]-lambda(3))*2*pi/lambda(2),(lambda(5)+1)/2)+lambda(4);
c=sumsqr((dat(:)-swth(:)));