% 
% [K,R] = gl2kr(A,G,C,L0)
% 
% Description:
%          Solve for the Kalman gain (K) and the innovation covariance (R)
%          The resulting model is of the form:
%           
%                 x_{k+1} = A x_k + K e_k
%                   y_k   = C x_k + e_k
%                cov(e_k) = R
%                
% Copyright: 
%          Peter Van Overschee, December 1995
%          peter.vanoverschee@esat.kuleuven.ac.be
%
%

function [K,R] =  gl2kr(A,G,C,L0)

if isempty(G) || isempty(L0)
  K = [];
  R = [];
else
  % Solve the Riccati equation
  [P,flag] = solvric(A,G,C,L0);
  if (flag == 1)
    disp('Warning: Non positive real covariance model => K = R = []');
    K = [];
    R = [];
  else
    % Make output (Page 63 for instance)
    R = L0 - C*P*C';
    K = (G - A*P*C')*inv(R);
  end
end


