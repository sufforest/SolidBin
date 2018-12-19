function [center, cluster_assignment, distortion]=kmeans_ortho_multiple(xx,k,num_ortho,num_random,threshold,max_iter)
% function [center, cluster_assignment, distortion]=kmeans_ortho_multiple(xx,k,num_ortho,num_random)
%
% This clusters the given n p-dim points xx(p,n) into k clusters. The
% algo makes num_ortho attempts by starting with orthogonal centers and
% then num_random attempts using random centers. It returns the centers
% and assignments with the minimum distortion. Each of  kmeans stops after the
% change is centers (norm of diff) is less than the 'threshold' or max_iterations.
  
  p=size(xx,1); 
  n=size(xx,2); 
  min_distortion=inf;
  total_num=num_random+num_ortho;
  
  norm_xx=normalize_2nd(xx); 
  for i=1:total_num

	%assign the center index ortho or randomly
	if (i <= num_ortho)
	  center_index=gen_orthogonal_centers(norm_xx); 
	else
	  is_ran_center=zeros(1,n+1); 
	  is_ran_center(n+1)=1;
	  center_index=zeros(1,k); 
	  for j=1:k
		ran_index=n+1;
		while is_ran_center(ran_index)
		  ran_index=floor(1+n*rand(1));
		end
		center_index(j)=ran_index;
	  end
	end  
	% run the kmeans on them
	kcenter=xx(:,center_index);
	[cluster_centers assignment single_distortion]=kmeans_nostat(max_iter,kcenter, threshold,xx);

	% see if distortion is reduced. 
	if single_distortion < min_distortion
	  min_distortion=single_distortion;
	  min_centers=cluster_centers; 
	  min_assignment=assignment;
	end
	
	center=min_centers; 
	cluster_assignment=min_assignment;
	distortion=min_distortion; 
  end
  
  
