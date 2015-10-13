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
end

idx=find((sign(x(1:end-1)).*sign(x(2:end)))<0)+1;

idx=[1; idx; numel(x)];

val=x(idx);
