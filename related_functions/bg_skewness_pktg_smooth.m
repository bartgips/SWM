function [skwIdx, brdOut, xOut]=bg_skewness_pktg_smooth(x)
% [skwIdx, brdOut, xOut]=bg_skewness_pktg_smooth(x)
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



% demean
x=bsxfun(@minus,x,nanmean(x));

%estimate base frequency
L=size(x,1);
nfft=2^nextpow2(L*10);
ftx=fft(x,nfft);
ftxpow=abs(ftx(1:(nfft/2+1)));
[~,periodEst]=max(ftxpow);
periodEst=nfft./periodEst;
skwIdx=nan(size(x,2),1);
brdOut=nan(size(x,2),3);

interpFac=100; %increase resolution 100-fold;

if nargout>2
  xOut=nan(L*interpFac,size(x,2));
end


for n=1:size(x,2)
  % interpolating
  t=1:L;
  tint=linspace(1,L,ceil(interpFac*L));
  xint=spline(t',x(:,n),tint');
  periodEst=periodEst*interpFac;
  dx=diff(xint);
  
  zeroCross=bg_find_zerocross(dx);
  pktgIdx=sign(diff(xint(zeroCross)));
  pktgIdx=[-pktgIdx(1,:); pktgIdx];
  
  % find the two most central extrema
  slopePos=1./abs(zeroCross-zeroCross(end)/2);
  [~,mLenidx]=max(conv(slopePos,ones(1,3),'same'));
  brd=zeroCross(mLenidx+[0:2]-1);
  sigma=0;
  breakloop=false;
  
  while diff(brd([1 3]))<.8*periodEst || breakloop
    
    sigma=sigma+interpFac/2;
    tau=-4*sigma:4*sigma;
    h=exp(-tau'.^2/(2*sigma^2));
    
    xsmth=conv2(xint,h,'same');
    dxsmth=diff(xsmth);
    dxsmth=dxsmth(1:end-1,:);
    
    zeroCross=bg_find_zerocross(dxsmth);
    % pktgIdx=sign(diff(xsmth(zeroCross)));
    % pktgIdx=[-pktgIdx(1,:); pktgIdx];
    
    % find the two most central extrema
    slopePos=1./abs(zeroCross-zeroCross(end)/2);
    [~,mLenidx]=max(conv(slopePos,ones(1,3),'same'));
    brd=zeroCross(mLenidx+[0:2]-1);
    
    if sigma>periodEst/2
      breakloop=true;
    end
    
%     figure(5)
%     plot(xsmth)
%     vline(brd)
%     drawnow
%     pause
  end
  xOut(:,n)=xsmth;
  
  slopeLen=diff(brd);
  slopeSign=sign(diff(xsmth(brd(1:2))));
  
  
  if slopeSign>0
    skwIdx(n)=(slopeLen(1)/sum(slopeLen)-.5)*2;
  elseif slopeSign<0
    skwIdx(n)=(slopeLen(2)/sum(slopeLen)-.5)*2;
  end
  
  brdOut(n,:)=brd;
end


