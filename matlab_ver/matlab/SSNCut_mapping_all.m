function uu=SSNCut_mapping(S,k,ML,CL,alpha,beta)
%function uu=SSNCut_mapping(S,k,ML,CL,alpha,beta)
%Returns the vectors using the SSNCut algorithm
%Input:
% S: similarity matrix
% k: number of clusters in the dataset
% ML:must-link constraint matrix
% CL:cannot-link constraint matrix
%alpha: parameter for must-link
%beta: parameter for cannot-link
%Output:
% uu: top k eigenvectors of the objective function
 
D=diag(sum(S)); 
Dsqrt=sqrt(D); 
W=S;                                        %original similarity matrix

if ~isempty(ML)&(alpha~=0)                  
   alpha=alpha*eye(size(S,1));
%   n_point = size(S,1);
%   num_ML= size(ML,1);
%   if(num_ML>15000)
%   ML0= zeros(n_point);
%   for i=1:n_point
 %      ML0(i,:)=ML(:,i)'*ML;
 %  end
  % W=W-alpha*ML'*ML;
%   else                          %add must-link constraint 
%       ML0=ML'*ML;
%   end
%   W= W-alpha*ML0;
W= W-alpha*(diag(sum(ML))-ML);%%%%%

end
        
if ~isempty(CL)&CL'==CL 
   if size(CL,1)~=1
       W=W-beta*CL;                           %add cannot-link constraint                      
   end
end

%if CL==CL' disp('CL ¶Ô³Æ');   end

L=Dsqrt\W/Dsqrt;
L= (L+L')/2;
%if L==L' disp('L ¶Ô³Æ');   end
%disp(sum(sum(L-L')));


[uu, dummy]=myeigs(L,k);
uu=uu'; %k*n

%uu=normalize_2nd(uu);

