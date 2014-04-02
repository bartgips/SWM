function [idx, val]=bg_find_zerocross(x)
% [idx, val]=bg_find_zerocross(x)
% 
% finds zerocrossings in x; useful when x is a derivative (to find local
% extrema)
% 
% Note: It always returns the first value after the zero-crossing, so the
% real 0 crossing happened somewhere between idx-1 and idx.
% 
% output:
% idx:    linear index of the zero crossing
% val:    the value of x closest to the zero crossing

if ~isvector(x)
  error('input should be a vector')
elseif isrow(x)
  x=x.';
end

valSign=sign(x(1));
idx=[];
lastZero=0;
done=0;
while ~done
  idxdum=[find(sign(x(lastZero+1:end))==-valSign,1)];
  if isempty(idxdum)
    done=1;
  else  
    idx=[idx idxdum+lastZero];
    lastZero=idxdum+lastZero;
    valSign=-valSign;
  end
end

idx=[1 idx numel(x)];

val=x(idx);
