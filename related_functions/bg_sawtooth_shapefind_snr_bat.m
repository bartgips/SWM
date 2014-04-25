function [skw]=bg_sawtooth_shapefind_snr_bat(snr,skew,tLen,numIt)
% [skw, stats]=bg_sawtooth_shapefind_snr_bat(snr,skew,tLen,numIt)

skw=nan(numel(snr), numel(skew), numel(tLen), numIt, 2);




for snrIt=1:numel(snr);
  for skwIt=1:numel(skew)
    for tLenIt=1:numel(tLen)
      for it=1:numIt
        skw(snrIt,skwIt,tLenIt,it,:)=bg_sawtooth_shapefind_snr(...
          snr(snrIt),skew(skwIt),tLen(tLenIt),0);
%         disp(it)
      end
    end
  end
end