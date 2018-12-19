function [cluster_centers,cluster_assignments,distortion]=kmeans_nostat(max_iter,kcenter, threshold,uu)
% function [cluster_centers,cluster_assignments,distortion]=kmeans_nostat(max_iter,kcenter, threshold,uu)
%
% Performs the K means on the data given the points uu as column
% vectors, kcenter as the initial centers. max_iter as the maximum number
% of iterations. The kmeans stops after the change is centers (norm of
% diff0) is less than the 'threshold'

%  K-means

  k=size(kcenter,2); %number of clusters. 
  n=size(uu,2) ; % number of points 
  kdist = zeros( k, n ); % distances of points from the clusters. 

  converged=0; 
  for iter = 1:max_iter;
	oldcenter=kcenter;
	d=0; 
	% compute distances to centers
	for ik = 1:k;
	  kdist( ik, : ) = sqrt( sum(( uu - repmat( kcenter(:,ik ), [1,n] )).^2));
	end;

	% compute the new cluster assignments 
	[ ddummy cluster_assignments ] = min( kdist );
	
	% recompute centers
	for ik = 1:k;
	  iik = find( cluster_assignments == ik );
	  if( length( iik) > 0 )
        kcenter( :, ik ) = sum( uu( :, iik ), 2 )/ length( iik );
	  end;
	end;
	% compare old center to new ones and stop in case similar
	diff0=norm(kcenter-oldcenter); 
       %disp(diff0);
	if diff0 < threshold
	  converged=1;
	  str=sprintf('Done with K-means. diff0=%.8g iter=%d',diff0,iter);
	  break; 
	end;
  end
  if (~converged)
	warning('K means did not coverge');
	str=sprintf('Done with K-means. diff0=%.8g iter=%d',diff0,iter);
	disp(str);
  end

  cluster_centers=kcenter;
  distortion=measure_distortion(uu,cluster_assignments); 
