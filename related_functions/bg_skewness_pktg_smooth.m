function [skwIdx, brdOut, xOut]=bg_skewness_pktg_smooth(x, bias, interpFac)
% [skwIdx, brdOut, xOut]=bg_skewness_pktg_smooth(x,bias)
% Calculates skewness by finding the extrema of a shape. This shape is
% first interpolated (cubic spline) 100-fold.
% The extrema are found by detecting zero-crossings in the derivative.
% The algorithm detects the dominant frequency first and then applies
% gaussian smoothing to remove higher-frequency noise, until the extrema
% spacing agrees with the dominant frequency.
% 
% %%%%%%%
% INPUT:
% %%%%%%%
% x:    Matrix containing the shapes of which the skewness should be
%       quantified. Matrix should be NxM where N is the number of
%       timepoints and M is the number of shapes.
% bias: (optional) Index biasing the algorithm to picking a period closest
%       to this index. (When the shape contains more than 1.5 period it
%       could detect pk-tg-pk or tg-pk-tg, the bias can force the algorithm
%       to always pick one of the two.)
%       (default= centre of window, i.e. N/2+.5);
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
brdOut=nan(size(x,2),3);

if nargin<3
interpFac=100; %increase resolution 100-fold;
end

if numel(bias)==1;
  bias=[bias-1 bias bias+1];
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
  dx=diff(xint);
  
  zeroCross=bg_find_zerocross(dx);
  zeroCrossOrig=zeroCross;
  pktgIdx=sign(diff(xint(zeroCross)));
  pktgIdx=[pktgIdx(1:end-1)];
  
  sigma=0;
  breakloop=false;
  meanperiod=zeroCross(2:end-1);
  meanperiod=(mean(diff(meanperiod(pktgIdx<0)))+mean(diff(meanperiod(pktgIdx>0))))/2;
    
  xsmth=xint;
  while meanperiod<.8*periodEstn && ~breakloop
    
    sigma=sigma+interpFac/2;
    tau=-4*sigma:4*sigma;
    h=exp(-tau'.^2/(2*sigma^2));
    h=h/sum(h);
    
    xsmth=conv2(xint,h,'same');
    dxsmth=diff(xsmth);
    dxsmth=dxsmth(1:end-1,:);
    
    zeroCross=bg_find_zerocross(dxsmth);
    
    pktgIdx=[nan; sign(diff(xsmth(zeroCross)))];
    
    % only consider peaks and troughs that are near the global max/min
    tUs=interpFac*5:interpFac*10:numel(xsmth);
    xUs=interp1(1:numel(xsmth),xsmth,tUs);
    
    dxUs=diff(xUs);
    dxUs=dxUs(1:end-1);
    zCUs=bg_find_zerocross(dxUs);
    zCUs=zCUs(2:end-1);
    
    
    
    zCUs=tUs(zCUs);
    rmsel=zeros(size(zeroCross));
    for k=1:numel(zCUs)
      rmsel=rmsel+(abs(zeroCross-zCUs(k))<(numel(tint)/numel(zCUs)*.1));
    end
    rmsel=logical(rmsel);  
    
    if sum(rmsel)<3 % if pruning was too rigorous, undo it
      rmsel=true(size(rmsel));
    end
    
    pktgIdx=pktgIdx(rmsel);
    zeroCross=zeroCross(rmsel);
    
    meanperiod=(mean(diff(zeroCross(pktgIdx<0)))+mean(diff(zeroCross(pktgIdx>0))))/2;
    
    
    if sigma>periodEstn/2
      breakloop=true;
    end 
    
%     figure(100)
%     plot(xsmth)
%     vline(zeroCross)
%     drawnow
%     pause
  end
  
  % find the period that's closest to the bias
%   zeroCross=zeroCross(2:end-1);
  dist=nan(numel(zeroCross)-2,1);
  for k=1:numel(zeroCross)-2
    dist(k)=sum(abs(zeroCross([0:2]+k)-bias));
  end  
  [~,mLenidx]=min(dist);
  brd=zeroCross(mLenidx+[0:2]);
  brd=zeroCrossOrig(bg_findnearest(zeroCrossOrig,brd));
  
  xOut(:,n)=xint;
  
  slopeLen=diff(brd);
  slopeSign=sign(diff(xsmth(brd(1:2))));
  
  
  if slopeSign>0
    skwIdx(n)=(slopeLen(1)/sum(slopeLen)-.5)*2;
  elseif slopeSign<0
    skwIdx(n)=(slopeLen(2)/sum(slopeLen)-.5)*2;
  end
  
  brdOut(n,:)=brd;
end


