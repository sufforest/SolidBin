function vv=normalize_2nd(vv)
% function vv=normalize_2nd(vv)
% normalized the vector vv along the second dimension. i.e. 
% vv returns has unit norm COLUMNS. 
% (Note: If the column is zero it is untouched and NO error msg is
% printed )   
  n=size(vv,2); 
  ss = sqrt( sum( vv.^2, 1 ));   % normalize
  for i=1:n
	if ss(i)==0 
	  vv(:,i)=0;
	else
	  vv(:,i)=vv(:,i)/ss(i);
	end
  end

  
