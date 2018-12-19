function distortion=measure_distortion(xx,ca)
% function distortion=measure_distortion(xx,ca)
%
% Measure the distortion of the points to the cluster centroids  w.r.t
% cluster assignment (ca). Points xx are assumed to be COLUMN
% vectors. The distortion is the variance w.r.t to the centroid. 

  distortion=0; 
  k=max(ca); 
  for i=1:k
	ci=find(ca==i); 
	len_i=length(ci);
	if (len_i > 0 ) 
	  points_i=xx(:,ci);
	  mean_i=mean(points_i,2); 
	  for j=1:len_i
		diff_ji=points_i(:,j)-mean_i; 
		distortion = distortion + (sum(diff_ji .^2));
	  end
	end
  end
  
