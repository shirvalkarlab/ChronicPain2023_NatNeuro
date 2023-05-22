% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Shanechi Lab, University of Southern California, 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcAIC Calculates the AIC for the given observation
%   AIC is calculated with this formula:
%       AIC = 2*nParam + N*log(det(cov(e, 1))) + N*ny*(1+log(2*pi))
%   This assumes that error e has multivariate Gaussian distribution.
%   With this assumption, the formula is derived from the loglikelihood 
%   of e given the model. In other words: 
%       log( P(e|model) ) = -N/2*log(det(cov(e, 1))) - N/2*ny*(1+log(2*pi))
%   and AIC = 2*nParam - 2*log( P(e|model) )%   Inputs:
%
%   Inputs:
%     - (1) y: the observation matrix. Rows are samples, columns are
%              dimensions of the data
%     - (2) yhat: prediction of the model. Same size as y.
%     - (3) nParam: number of free parameters in the model
%   Note: For AIC from several models to be comparable, they must have all
%         used the exact same data y
%   Outputs:
%     - (1) AIC calculated based on Gaussian error assumption
%   Usage example:
%       % This example uses functions subid and predic from:
%       https://www.mathworks.com/matlabcentral/fileexchange/2290-subspace-identification-for-linear-systems
%       % Identify an LSSM for the data y
%       horiz = 10; % supspace system identification horizon
%       ord = 15; % Model order
%       [A,B,C,D,K,R,AUX] = subid(y,[],horiz,ord,[],[],1); % 
%       [yp,erp] = predic(y,[],A,B,C,D,K); % Run kalman over data using this LSSM
%       nParam = (2*ny+1)*ord; % Number of estimated parameters in LSSM
%       AIC = calcAIC(y, yp, nParam);
  
function AIC = calcAIC(y, yhat, nParam)

if any(size(yhat)~=size(y)), error('Dimension mismatch!'); end
if (size(y, 1) < size(y, 2)), warning('Not enough observations. Error Covariance will be singular -> AIC = -Inf' ); end
N = size(y, 1); % Num of observations
ny = size(y, 2); % Dimension of observations
e = y - yhat;
eCOV = cov(e, 1, 'omitrows');
% Estimate error covariance
AIC1 = 2*nParam + N*log(det(eCOV)) + N*ny*(1+log(2*pi)); % THIS CAN BE BEYOND NUMERICAL ACCURACY because of 'det' operation
if any(isnan(eCOV(:))) || any(isnan(e(:)))
    AIC = nan;
elseif any(isinf(eCOV(:)))
    AIC = -Inf;
else
    % logDetCOV = sum(log(eig(eCOV))); % not totally numerically stable
    % Since COV must be symmetric and positive semi-definite, we can use 
    % SVD to compute its singular values
    [~,S,~] = svd(eCOV,'econ');
    logDetCOV = sum(log(diag(S)));
    AIC = 2*nParam + N*logDetCOV + N*ny*(1+log(2*pi));
end

if (abs(AIC1)~=Inf) && abs(AIC1-AIC)>0.001*abs(AIC); fprintf('WARNING: Numerical accuracy of AIC computation need attention!\n'); end

end