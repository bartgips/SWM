function [stats]=bg_jackknife_sawtooth(shapeMat,verbose,fignum)
% stats = bg_jackknife_sawtooth(shapeMat, verbose, fignum)
%
% Estimates confidence interval of skewness of noisy shapes contained in
% 'shapeMat' using delete-1 jackknife. This means the number of samples used
% for estimation is equal to the number of shapes in shapeMat.
%
% This method is useful when the individual shapes are very
% noisy, and this noise hinders you in determining the true skewness of
% every shape individually. Averaging over all-but-one shapes, will average
% out the noise as much as possible.
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

if nargin<2
  verbose=1;
end
reverseStr=[];

if nargin >2 && fignum
  figure(fignum)
  clf
  set(fignum,'visible','off')
end

[numTemp, tempLen]=size(shapeMat);
numIt=numTemp;
skwIdx=nan(numIt,1);

% cut off sides (focus on centre theta shape +- 10%)
[~,ax]=findpeaks(sign(nansum(shapeMat(:)-mean(minmax(nanmean(shapeMat)))))*nanmean(shapeMat),'minpeakdistance',round(tempLen/3));
brd=mean(ax)+[-1 1]*diff(ax)*.55;
brd=[max(floor(brd(1)),1) min(ceil(brd(2)),tempLen)];
shapeMat=shapeMat(:,brd(1):brd(2));

% intial fit parameters
X0=[diff(minmax(nanmean(shapeMat)))/2, tempLen/2, 0, mean(shapeMat(:)), 0]';
fitOptions=optimset('Algorithm','interior-point','Display','off');

for iter=1:numIt
  
  sel=true(size(shapeMat,1),1);
  sel(iter)=false;
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
  
  if nargin >2 && fignum
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
stats.sem=sqrt((numTemp-1)*var(skwIdx,1));
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
