function Predict_Labels = SSNCut(Simi_Matrix, K_Clust, ML, CL, alpha, beta)

%function SSNCut(Simi_Matrix, K_Clust)
%Cluster using the semi-supervised normalized cut method SSNCut 
%Input:  
% Simi_Matrix: similarity between docs of the dataset
% K_Clust: cluster number of the dataset
% ML: constraints of must-link pairs
% CL: constraints of cannot-link pairs
% alpha: paramater for positive constraint
% beta: parameter for negative constraint
%Output:
%predict_Labels: result Labels for docs in the dataset

global_options;
% compute the top k eigenvectors of the the cost function of SSNCut  
xx=SSNCut_mapping_all(Simi_Matrix,K_Clust, ML, CL, alpha, beta);  

global KMEANS_THRESHOLD KMEANS_MAX_ITER;
KMEANS_MAX_ITER=200;
threshold = KMEANS_THRESHOLD;
max_iter  = KMEANS_MAX_ITER;

%Clusters the given POINTS using kmeans
Predict_Labels=cluster_point_kmeans(xx,K_Clust,5,20); 
