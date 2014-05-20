function [skw, meanShape]=bg_sawtooth_shapefind_snr_bat(snr,skew,tLen,numIt,verbose)
% [skw, meanShape]=bg_sawtooth_shapefind_snr_bat(snr,skew,tLen,numIt)

skw=nan(numel(snr), numel(skew), numel(tLen), numIt, 2);

meanShape=nan(numel(snr), numel(skew), numel(tLen), numIt, 1000, 2);

if nargin<5
  verbose=false;
end

for snrIt=1:numel(snr);
  for skwIt=1:numel(skew)
    for tLenIt=1:numel(tLen)
      for it=1:numIt
        [skw(snrIt,skwIt,tLenIt,it,:), stats]=bg_sawtooth_shapefind_snr(...
          snr(snrIt),skew(skwIt),tLen(tLenIt),verbose);
        meanShape(snrIt,skwIt,tLenIt,it,:,:)=stats.meanShape;
      end
    end
  end
end