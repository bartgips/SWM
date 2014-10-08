function [skwIdx, brdOut, xOut]=bg_skewness_pktg_harsh(x, bias, interpFac, numExtr)
% [skwIdx, brdOut, xOut]=bg_skewness_pktg_harsh(x, bias, interpFac, numExtr)
% Calculates skewness by finding the extrema of a shape. This shape is
% first interpolated (cubic spline) with factor "interpFac".
% The extrema are found by detecting zero-crossings in the derivative.
% The algorithm detects the dominant frequency first and then applies
% gaussian smoothing to remove higher-frequency noise, until the extrema
% spacing agrees with the dominant frequency.
%
% %%%%%%%
% INPUT:
% %%%%%%%
% x:        Matrix containing the shapes of which the skewness should be
%           quantified. Matrix should be NxM where N is the number of
%           timepoints and M is the number of shapes.
% bias:     (optional) Index biasing the algorithm to picking a period closest
%           to this index. (When the shape contains more than 1.5 period it
%           could detect pk-tg-pk or tg-pk-tg, the bias can force the algorithm
%           to always pick one of the two.)
%           (default= centre of window, i.e. N/2+.5);
% interpFac:Interpolation factor. How many times the resolution in x will
%           be enhanced
% numExtr:  Number of extrema given in brdOut.
%           Default and minimum value is 3.
%
% %%%%%%%
% OUTPUT:
% %%%%%%%
% skwIdx: The skewness index for every of the M shapes.
%         The skewness index is a number between -1 and 1. -1 corresponds
%         to a sawtooth leaning completely to the left (i.e. the upgoing
%         flank is a heavside like function). +1 corresponds to a sawtooth
%         leaning to the left.
%
% brdOut: The detected 3 extrema to calculate the skewness index
% xOut:   (Optional) output of the interpolated and smoothed version of x
%         used for dermining the skewness.
if isrow(x)
  x=x(:);
end

if nargin<2 || isempty(bias)
  bias=size(x,1)/2+.5;
end

if nargin < 4
  numExtr=3;
end

% demean
x=bsxfun(@minus,x,nanmean(x));

%estimate base frequency
L=size(x,1);
nfft=2^nextpow2(L*10);
ftx=fft(x,nfft);
ftxpow=abs(ftx(1:(nfft/2+1),:));
[~,periodEst]=max(ftxpow);
periodEst=nfft./periodEst;
skwIdx=nan(size(x,2),1);
brdOut=nan(size(x,2),numExtr);

if nargin<3
  interpFac=100; %increase resolution 100-fold;
end

if numel(bias)==1;
  bias=[1:numExtr]+bias-(numExtr+1)/2;
end

if numel(bias) ~= numExtr
  bias=[1:numExtr]+mean(bias)-(numExtr+1)/2;
end

bias=bias*interpFac;

if nargout>2
  xOut=nan(L*interpFac,size(x,2));
end


for n=1:size(x,2)
  % interpolating
  t=1:L;
  tint=linspace(1,L,ceil(interpFac*L));
  xint=spline(t',x(:,n),tint');
  periodEstn=periodEst(n)*interpFac;
  [zeroCross, pktgIdx]=harshReduction(xint);
  
  sigma=0;
  breakloop=false;
  meanperiod=zeroCross;
  meanperiod=(mean(diff(meanperiod(pktgIdx<0)))+mean(diff(meanperiod(pktgIdx>0))))/2;
  
  xsmth=xint;
  while meanperiod<.8*periodEstn && ~breakloop
    
    [zeroCross, pktgIdx]=harshReduction(xsmth);

    meanperiod=(mean(diff(zeroCross(pktgIdx<0)))+mean(diff(zeroCross(pktgIdx>0))))/2;
    
    
    if sigma>periodEstn/2
      breakloop=true;
    end
    
    %     figure(100)
    %     plot(xsmth)
    %     vline(zeroCross)
    %     drawnow
    %     pause
    
    % increase smoothing.
    sigma=sigma+interpFac/2;
    tau=-4*sigma:4*sigma;
    h=exp(-tau'.^2/(2*sigma^2));
    h=h/sum(h);
    
    xsmth=conv2(xint,h,'same');
  end
  
  % find the period that's closest to the bias
  %   zeroCross=zeroCross(2:end-1);
  brd=nan(numExtr,1);
  if numel(zeroCross)<numExtr
    numExtr=numel(zeroCross);
    bias=bias(1:numExtr);
  end
  
  if numel(zeroCross)<3
    error('not enough extrema to calculate skewness')
  end
  
  
  
  dist=nan(numel(zeroCross)-2,1);
  for k=1:numel(zeroCross)-numExtr+1
    dist(k)=sum(abs(zeroCross([0:numExtr-1]+k)-bias));
  end
  [~,mLenidx]=min(dist);
  brd(1:numExtr)=zeroCross(mLenidx+[0:numExtr-1]);
  %   brd=zeroCrossOrig(bg_findnearest(zeroCrossOrig,brd));
  xOut(:,n)=xint;
  
  slopeLen=diff(brd);
  slopeSign=sign(diff(xsmth(brd(1:2))));
  
  
  if slopeSign>0 % first slope is rising
    t_rise=nanmean(slopeLen(1:2:end));
    t_fall=nanmean(slopeLen(2:2:end));
  elseif slopeSign<0 % first slope is descending
    t_rise=nanmean(slopeLen(2:2:end));
    t_fall=nanmean(slopeLen(1:2:end));
  end
  
  skwIdx(n)=(t_rise/(t_rise+t_fall) -.5) *2;
  
  brdOut(n,:)=brd;
end
end

function [zeroCross, pktgIdx]=harshReduction(xint)
dx=diff(xint);

zeroCross=bg_find_zerocross(dx);
pktgIdx=sign(diff(xint(zeroCross)));
pktgIdx=[pktgIdx(1:end-1)];
zeroCross=zeroCross(2:end-1);

% firstly, remove local extrema that are obvious
% (peaks that are lower than median trough and vice versa)
badSel=1;
while sum(badSel)>0
  pkQ=median(xint(zeroCross(pktgIdx==1)));
  tgQ=median(xint(zeroCross(pktgIdx==-1)));
  badPk=(xint(zeroCross)<tgQ) & pktgIdx==1;
  badTg=xint(zeroCross)>pkQ & pktgIdx==-1;
  badSel=badPk | badTg;
  zeroCross=zeroCross(~badSel);
  pktgIdx=pktgIdx(~badSel);
end

% Next pass, more agressive, use max and min instead of median
badSel=1;
while sum(badSel)>0
  pkQ=quantile(xint(zeroCross(pktgIdx==1)),0.0);
  tgQ=quantile(xint(zeroCross(pktgIdx==-1)),1.0);
  badPk=(xint(zeroCross)<tgQ) & pktgIdx==1;
  badTg=xint(zeroCross)>pkQ & pktgIdx==-1;
  if diff(minmax(xint(zeroCross(pktgIdx==1))))> 2* diff(minmax(xint(zeroCross(pktgIdx==-1))))
    badSel=badPk; % only remove bad peaks first, since the sprea in peak positions is too large
  elseif diff(minmax(xint(zeroCross(pktgIdx==1)))) < 1/2*diff(minmax(xint(zeroCross(pktgIdx==-1))))
    badSel=badTg;
  else
    badSel=badPk | badTg; % remove bad peaks and troughs at the same time
  end
  zeroCross=zeroCross(~badSel);
  pktgIdx=pktgIdx(~badSel);
end

% Finally remove troughs that are closer to median of peaks than median of
% troughs and  vice versa
badSel=1;
while sum(badSel)>0
  pkQ=quantile(xint(zeroCross(pktgIdx==1)),0.5);
  tgQ=quantile(xint(zeroCross(pktgIdx==-1)),0.5);
  badPk=abs(xint(zeroCross)-tgQ)<abs(xint(zeroCross)-pkQ) & pktgIdx==1;
  badTg=abs(xint(zeroCross)-tgQ)>abs(xint(zeroCross)-pkQ) & pktgIdx==-1;
  badSel=badPk | badTg;
  zeroCross=zeroCross(~badSel);
  pktgIdx=pktgIdx(~badSel);
end

% now remove consecutive troughs or peaks
breakLoop=false;
pktgIdx=[pktgIdx; nan]; % add NaN as a marker for the end of the vector
zeroCross=[zeroCross nan];
while ~breakLoop
  firstTg=find(pktgIdx==-1,1);
  firstPk=find(pktgIdx==1,1);
  firstNan=find(isnan(pktgIdx),1);
  
  %       figure(1); clf; plot(xint); hold on; plot((zeroCross(pktgIdx==1)),xint(zeroCross(pktgIdx==1)),'or'); plot((zeroCross(pktgIdx==-1)),xint(zeroCross(pktgIdx==-1)),'xk')
  %       vline((zeroCross([firstTg firstPk])))
  %       pause
  
  selVec=false(size(pktgIdx));
  if firstNan<firstTg && firstNan<firstPk % loop is done
    zeroCross=zeroCross(~isnan(pktgIdx));
    pktgIdx=pktgIdx(~isnan(pktgIdx));
    breakLoop=true;
  elseif firstNan<firstPk % loop almost done, final bit is all troughs
    selVec(firstTg:firstNan-1)=true;
    Seq=zeroCross(selVec);
    Seqpktg=pktgIdx(selVec);
    [~,mIdx]=min(xint(zeroCross(selVec)));
    zeroCross=[zeroCross(~selVec) Seq(mIdx)];
    pktgIdx=[pktgIdx(~selVec); Seqpktg(mIdx)];
  elseif firstNan<firstTg % loop almost done, final bit is all peaks
    selVec(firstPk:firstNan-1)=true;
    Seq=zeroCross(selVec);
    Seqpktg=pktgIdx(selVec);
    [~,mIdx]=max(xint(zeroCross(selVec)));
    zeroCross=[zeroCross(~selVec) Seq(mIdx)];
    pktgIdx=[pktgIdx(~selVec); Seqpktg(mIdx)];
  elseif (firstPk-firstTg)>1 %first extrema are troughs
    selVec(firstTg:firstPk-1)=true;
    Seq=zeroCross(selVec);
    Seqpktg=pktgIdx(selVec);
    [~,mIdx]=min(xint(zeroCross(selVec)));
    zeroCross=[zeroCross(~selVec) Seq(mIdx)];
    pktgIdx=[pktgIdx(~selVec); Seqpktg(mIdx)];
  elseif (firstTg-firstPk)>1 %first extrema are peaks
    selVec(firstPk:firstTg-1)=true;
    Seq=zeroCross(selVec);
    Seqpktg=pktgIdx(selVec);
    [~,mIdx]=max(xint(zeroCross(selVec)));
    zeroCross=[zeroCross(~selVec) Seq(mIdx)];
    pktgIdx=[pktgIdx(~selVec); Seqpktg(mIdx)];
  elseif abs(firstPk-firstTg)==1 % a single extremum: skip it
    zeroCross=zeroCross([2:end 1]);
    pktgIdx=pktgIdx([2:end 1]);
  end
end
end

