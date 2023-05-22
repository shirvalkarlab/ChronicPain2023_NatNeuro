% 
% [P,flag] = solvric(A,G,C,L0)
% 
% Description:
%       Solves the Forward Riccati equation:
%          
%       P = A P A' + (G - A P C') (L0 - C P C')^{-1} (G - A P C')' 
%       
%       Using the generalized eigenvalue decomposition of page 62       
%       flag = 1 when the covariance sequence is not positive real
%
%       To solve the backward Riccati equation:
%          
%       N = A' N A + (C' - A' N G) (L0 - G' N G)^{-1} (C' - A' N G)'
%       
%       Use:
%       [N,flag] = solvric(A',C',G',L0')
%                
%
% Copyright: 
%          Peter Van Overschee, December 1995
%          peter.vanoverschee@esat.kuleuven.ac.be
%
%

function [P,flag] = solvric(A,G,C,L0)

if isempty(L0) | isempty(G);
  P = [];
  flag = 0;
else
  [n,n]=size(A); 			% Dimensions
  L0i=inv(L0); 				% Compute the inverse once

  % Set up the matrices for the eigenvalue decomposition
  AA = [A'-C'*L0i*G' zeros(n,n);-G*L0i*G' eye(n)];
  BB = [eye(n) -C'*L0i*C;zeros(n,n) A-G*L0i*C];

  % Compute the eigenvalue decomposition
  [v,d] = eig(AA,BB);
  ew = diag(d);

  % If there's an eigenvalue on the unit circle => no solution
  flag = 0;
  if all(abs(abs(ew)-ones(2*n,1)) > 1e-9) < 1;flag = 1;end

  % Sort the eigenvalues
  [dum,I]=sort(abs(ew));

  % Compute P
  P=real(v(n+1:2*n,I(1:n))/v(1:n,I(1:n)));
end


