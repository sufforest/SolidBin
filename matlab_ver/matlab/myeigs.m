function [v,d]=myeigs(S,k) 
%
% function [v,d]=myeigs(S,k) 
% Note: d returned is a vector not a diag matrix 
% and may not be sorted ...
  opts.disp=0; 
  Sdiff=S-S'; 
  if max(max(abs(Sdiff))) < 1e-20 % the matrix is symm 
	[v d]=eigs(S,k,'la',opts); 
  else
	[v d]=eigs(S,k,'lr',opts); 
  end
  d=diag(d); 
  
	
  
