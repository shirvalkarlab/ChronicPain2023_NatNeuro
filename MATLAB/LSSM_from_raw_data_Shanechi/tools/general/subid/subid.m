% 
%   General subspace identification 
%   -------------------------------
%   
%   The algorithm 'subid' identifies deterministic, stochastic 
%   as well as combined state space systems from IO data.
%
%           [A,B,C,D,K,R] = subid(y,u,i);
% 
%   Inputs:
%           y: matrix of measured outputs;
%                   - First dim must be time, and 2nd dim output dimensions
%           u: matrix of measured inputs;
%                   - First dim must be time, and 2nd dim output dimensions
%              for stochastic systems, u = []
%           i: number of block rows in Hankel matrices
%              (i * #outputs) is the max. order that can be estimated 
%              Typically: i = 2 * (max order)/(#outputs)
%           
%   Outputs:
%           A,B,C,D,K,R: combined state space system
%           
%                  x_{k+1) = A x_k + B u_k + K e_k        
%                    y_k   = C x_k + D u_k + e_k
%                 cov(e_k) = R
%                 
%           For deterministic systems: K = R = []
%           For stochastic systems:    B = D = []
%
%   Optional:
%
%           [A,B,C,D,K,R,AUX,ss] = subid(y,u,i,n,AUX,W,sil);
%   
%           n:    optional order estimate (default [])
%                 if not given, the user is prompted for the order
%           AUX:  optional auxilary variable to increase speed (default [])
%           W:    optional weighting flag
%                      SV:    Singular values based algorithm 
%                             (default for systems with input u)
%                      CVA:   Canonical variate based algorithm
%                             (default for systems without input u)
%           ss:   column vector with singular values
%           sil:  when equal to 1 no text output is generated
%           forceStable (default: false): if true will chnage steps 4 and 5
%           according to book (page 129) to ensure stability.
%           
%   Example:
%   
%           [A,B,C,D,K,R,AUX] = subid(y,u,10,2);
%           for k=3:6
%              [A,B,C,D] = subid(y,u,10,k,AUX);
%           end
%           
%   Reference:
%   
%           Subspace Identification for Linear Systems
%           Theory - Implementation - Applications
%           Peter Van Overschee / Bart De Moor
%           Kluwer Academic Publishers, 1996
%           Stochastic algorithm:   Figure 3.13 page 90 (positive)
%           Combined algorithm:     Figure 4.8 page 131 (robust)
%
%   Copyright:
%   
%           Peter Van Overschee, December 1995
%           peter.vanoverschee@esat.kuleuven.ac.be
%
%

function [A,B,C,D,K,Ro,AUX,ss,Qs,Rs,Ss] = subid(y,u,i,n,AUXin,W,sil,forceStable)

if nargin < 8, forceStable = false; end

warning off
  
if (nargin < 7);sil = 0; end

mydisp(sil,' ');
mydisp(sil,'   Subspace Identification');
mydisp(sil,'   -----------------------');

% Check the arguments
if (nargin < 3);error('subid needs at least three arguments');end
if (nargin < 4);n = [];end
if (nargin < 5);AUXin = [];end

% Check if its deterministic or stochastic ID
if isempty(u);   ds_flag = 2; 		% Stochastic
else;           ds_flag = 1; 		% Deterministic
end  

% Give W its default value
if (nargin < 6);W = [];end
if isempty(W)
  if (ds_flag == 1); W = 'SV'; 		% Deterministic: default to SV
  else;            W = 'CVA';end 	% Stochastic: default to CVA
end


% Turn the data into row vectors and check

[l,ny] = size(y);if (ny < l);y = y';[l,ny] = size(y);end
segCount = 1;


if (i < 0);error('Number of block rows should be positive');end
if (l < 0);error('Need a non-empty output vector');end
if (ds_flag == 1)
    
  [m,nu] = size(u);if (nu < m);u = u';[m,nu] = size(u);end
  
  if (m < 0);error('Need a non-empty input vector');end
  if (nu ~= ny);error('Number of data points different in input and output');end
else
  m = 0;
end

% Determine the number of columns in the Hankel matrices
if segCount == 1
    j = ny-2*i+1;
else
    j = ny-(2*i-1)*segCount;
end

if (j < (2*l*i));error('Not enough data points');end

% Check the weight to be used
Wn = 0;
if (length(W) == 2) 
  if (all(W == 'SV') | all(W == 'sv') | all(W == 'Sv'));
    Wn = 1; 
    if (ds_flag == 1);Waux = 2;else;Waux = 3;end
  end
end    
if (length(W) == 3) 
  if (prod(W == 'CVA') | prod(W == 'cva') | prod(W == 'Cva'));
    Wn = 2;
    if (ds_flag == 1);Waux = 3;else;Waux = 1;end
  end 
end
if (Wn == 0);error('W should be SV or CVA');end
W = Wn;

% Check compatibility of AUXin
if (ds_flag == 1)
  if ~iscell(u)
    Uaux = u(1,1);
  else
    Uaux = u{1}(1,1);
  end
else;Uaux = [];end
if ~iscell(y)
  Yaux = y(1,1);
else
  Yaux = y{1}(1,1);
end
[AUXin,Wflag] = chkaux(AUXin,i,Uaux,Yaux,ds_flag,Waux,sil);

  
% Compute the R factor
if isempty(AUXin)
  if ~iscell(y)
    Y = blkhank(y/sqrt(j),2*i,j); 	% Output block Hankel
  else
    Y = [];
    for si = 1:length(y)
      thisJ = size(y{si}, 2)-2*i+1;
      Y = cat(2, Y, blkhank(y{si}/sqrt(j),2*i,thisJ,false)); 	% Output block Hankel
    end
  end
  mydisp(sil,'      Computing ... R factor');
  if (ds_flag == 1)
    if ~iscell(u)
      U = blkhank(u/sqrt(j),2*i,j); 	% Input block Hankel
    else
      U = [];
      for si = 1:length(u)
        thisJ = size(u{si}, 2)-2*i+1;
        U = cat(2, U, blkhank(u{si}/sqrt(j),2*i,thisJ)); 	% Input block Hankel
      end
    end
    R = triu(qr([U;Y]'))'; 		% R factor
    clear U Y
  else
    R = triu(qr(Y'))'; 			% R factor
    clear Y
  end
  R = R(1:2*i*(m+l),1:2*i*(m+l)); 	% Truncate
else
  R = AUXin(2:2*i*(m+l)+1,1:2*(m+l)*i);
  bb = 2*i*(m+l)+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              BEGIN ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% **************************************
%               STEP 1 
% **************************************

mi2  = 2*m*i;
% Set up some matrices
if isempty(AUXin) | (Wflag == 1)
  Rf = R((2*m+l)*i+1:2*(m+l)*i,:); 	% Future outputs
  Rp = [R(1:m*i,:);R(2*m*i+1:(2*m+l)*i,:)]; % Past (inputs and) outputs
  if (ds_flag == 1)
    Ru  = R(m*i+1:2*m*i,1:mi2); 	% Future inputs
    % Perpendicular Future outputs 
    Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)]; 
    % Perpendicular Past
    Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)]; 
  end
end

% The oblique projection:
% Computed as 6.1 on page 166
% obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
% The extra projection on Ufp (Uf perpendicular) tends to give better
% numerical conditioning (see algo page 131)
% And it is needed both for CVA as MOESP

if isempty(AUXin)
  if (ds_flag == 1)
    % Funny rank check (SVD takes too long)
    % This check is needed to avoid rank deficiency warnings
    if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10
      Ob  = (Rfp*pinv(Rpp')')*Rp; 	% Oblique projection
    else
      Ob = (Rfp/Rpp)*Rp;
    end
  else    
    % Ob  = (Rf/Rp)*Rp; which is the same as 
    Ob = [Rf(:,1:l*i),zeros(l*i,l*i)];
  end
else
  % Determine Ob from AUXin
  Ob = AUXin(bb+1:bb+l*i,1:2*(l+m)*i);
  bb = bb+l*i;
end


% **************************************
%               STEP 2 
% **************************************

% Compute the SVD
if isempty(AUXin) | (Wflag == 1)
  mydisp(sil,'      Computing ... SVD');
  % Compute the matrix WOW we want to take an SVD of
  % W = 1 (SV), W = 2 (CVA)
  if (ds_flag == 1)
    % Extra projection of Ob on Uf perpendicular
    WOW = [Ob(:,1:mi2) - (Ob(:,1:mi2)/Ru)*Ru,Ob(:,mi2+1:2*(m+l)*i)];
  else
    WOW = Ob;
  end    
  if (W == 2)
    W1i = triu(qr(Rf'));
    W1i = W1i(1:l*i,1:l*i)';
    WOW = W1i\WOW;
  end
  [U,S,V] = svd(WOW);
  if W == 2;U = W1i*U;end 		% CVA
  ss = diag(S);
  clear V S WOW
else
  U = AUXin(bb+1:bb+l*i,1:l*i);
  ss = AUXin(bb+1:bb+l*i,l*i+1);
end


% **************************************
%               STEP 3 
% **************************************

% Determine the order from the singular values
if isempty(n)
  figure(gcf);hold off;subplot;
  if (W == 2)
    bar([1:l*i],real(acos(ss))*180/pi);
    title('Principal Angles');
    ylabel('degrees');
  else
    H = bar([1:l*i],ss); xx = get(H,'XData'); yy = get(H,'YData'); 
    semilogy(xx,yy+10^(floor(log10(min(ss)))));
    axis([0,length(ss)+1,10^(floor(log10(min(ss)))),10^(ceil(log10(max(ss))))]);
    title('Singular Values');
  end
  xlabel('Order');
  n = 0;
  while (n < 1) | (n > l*i-1)
    n = input('      System order ? ');
    if isempty(n);n = -1;end
  end
  mydisp(sil,' ')
end

U1 = U(:,1:n); 				% Determine U1


% **************************************
%               STEP 4 
% **************************************

% Determine gam and gamm
gam  = U1*diag(sqrt(ss(1:n)));
gamm = gam(1:l*(i-1),:);
% The pseudo inverses
gam_inv  = pinv(gam); 			% Pseudo inverse
gamm_inv = pinv(gamm); 			% Pseudo inverse

if forceStable % Modified according to pvo book page 129
    % Add zeros to gamm to get gam0 (will still call it gamm here)
    gamm0 = [gam(l+(1:l*(i-1)),:); zeros(l, size(gam, 2))];
end

% **************************************
%               STEP 5 
% **************************************


% Determine the matrices A and C
mydisp(sil,['      Computing ... System matrices A,C (Order ',num2str(n),')']); 
Rhs = [  gam_inv*R((2*m+l)*i+1:2*(m+l)*i,1:(2*m+l)*i),zeros(n,l) ; ...
    R(m*i+1:2*m*i,1:(2*m+l)*i+l)];
Lhs = [        gamm_inv*R((2*m+l)*i+l+1:2*(m+l)*i,1:(2*m+l)*i+l) ; ...
    R((2*m+l)*i+1:(2*m+l)*i+l,1:(2*m+l)*i+l)];

if ~forceStable % Original
    % Solve least square
    sol = Lhs/Rhs;

    % Extract the system matrices A and C
    A = sol(1:n,1:n);
    C = sol(n+1:n+l,1:n);
    res = Lhs - sol*Rhs; 			% Residuals
else % Modified according to pvo book page 129
    % Extract the system matrices A and C
    C = gam(1:l, :);
    A = gam_inv*gamm0;
    res = Lhs - [A;C]*Rhs; 			% Residuals
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Recompute gamma from A and C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam=C;
for k=2:i
	gam((k-1)*l+1:k*l,:) = gam((k-2)*l+1:(k-1)*l,:)*A;
end
gamm = gam(1:l*(i-1),:);      
gam_inv = pinv(gam);
gamm_inv = pinv(gamm);	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Recompute the states with the new gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rhs = [  gam_inv*R((2*m+l)*i+1:2*(m+l)*i,1:(2*m+l)*i),zeros(n,l) ; ...
    R(m*i+1:2*m*i,1:(2*m+l)*i+l)];
Lhs = [        gamm_inv*R((2*m+l)*i+l+1:2*(m+l)*i,1:(2*m+l)*i+l) ; ...
    R((2*m+l)*i+1:(2*m+l)*i+l,1:(2*m+l)*i+l)];


% **************************************
%               STEP 6 
% **************************************

if (ds_flag == 2)
  B = [];
  D = [];
else
  mydisp(sil,['      Computing ... System matrices B,D (Order ',num2str(n),')']); 
  % P and Q as on page 125
  P = Lhs - [A;C]*Rhs(1:n,:);
  P = P(:,1:2*m*i);
  Q = R(m*i+1:2*m*i,1:2*m*i); 		% Future inputs

  % L1, L2, M as on page 119
  L1 = A * gam_inv;
  L2 = C * gam_inv;
  M  = [zeros(n,l),gamm_inv];
  X  = [eye(l),zeros(l,n);zeros(l*(i-1),l),gamm];
  
  totm=0;
  for k=1:i
    % Calculate N and the Kronecker products (page 126)
    N = [...
	    [M(:,(k-1)*l+1:l*i)-L1(:,(k-1)*l+1:l*i),zeros(n,(k-1)*l)]
	[-L2(:,(k-1)*l+1:l*i),zeros(l,(k-1)*l)]];
    if k == 1;
      N(n+1:n+l,1:l) = eye(l) + N(n+1:n+l,1:l);
    end
    N = N*X;
    totm = totm + kron(Q((k-1)*m+1:k*m,:)',N);
  end
  
  % Solve Least Squares
  P = P(:);
  sol = totm\P;
  
  % Find B and D
  sol_bd = reshape(sol,(n+l),m);
  D = sol_bd(1:l,:);
  B = sol_bd(l+1:l+n,:);
end  


% **************************************
%               STEP 7 
% **************************************

if (norm(res) > 1e-10)
  % Determine QSR from the residuals
  mydisp(sil,['      Computing ... System matrices G,L0 (Order ',num2str(n),')']); 
  % Determine the residuals
  cov = res*res'; 			% Covariance
  Qs = cov(1:n,1:n);Ss = cov(1:n,n+1:n+l);Rs = cov(n+1:n+l,n+1:n+l); 
  
  sig = dlyap(A,Qs);
  G = A*sig*C' + Ss;
  L0 = C*sig*C' + Rs;

  % Determine K and Ro
  mydisp(sil,'      Computing ... Riccati solution')
  [K,Ro] = gl2kr(A,G,C,L0);
else
  Ro = [];
  K = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  END ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Make AUX when needed
if nargout > 6
  AUX = zeros((4*l+2*m)*i+1,2*(m+l)*i);
  if isempty(Uaux);Uaux = 0;end
  info = [ds_flag,i,Uaux,Yaux,Waux]; % in/out - i - u(1,1) - y(1,1) - W
  AUX(1,1:5) = info;
  bb = 1;
  AUX(bb+1:bb+2*(m+l)*i,1:2*(m+l)*i) = R;
  bb = bb+2*(m+l)*i;
  AUX(bb+1:bb+l*i,1:2*(m+l)*i) = Ob;
  bb = bb+l*i;
  AUX(bb+1:bb+l*i,1:l*i) = U;
  AUX(bb+1:bb+l*i,l*i+1) = ss;
end