function [val, xOut]=isfieldi(x,fieldName)
% checks wheter structure x contains field fieldName independent of the
% case (lower or upper). If x contains such a field, but the spelling is
% off, the output x will have it's field renamed

if iscell(fieldName)
  val=nan(size(fieldName));
  for n=1:numel(fieldName)
    % recursively call for every cell
    [val(n),x]=isfieldi(x,fieldName{n});
  end
  xOut=x;
else
  names=fieldnames(x);
  selVec=strcmpi(fieldName,names);
  if sum(selVec)>1
    selVec=strcmp(fieldName,names);
    if ~any(selVec)
      error('structure array contains multiple matches, but no exact one')
    else
      val=true;
      xOut=x;
    end
  elseif any(selVec)
    val=true;
    oldField=names{selVec};
    fieldContents=x.(oldField);
    xOut=rmfield(x,oldField);
    xOut.(fieldName)=fieldContents;
  else
    val=false;
    xOut=x;
  end
end

