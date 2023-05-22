% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Yuxiao Yang, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2018
%   Code from: https://www.nature.com/articles/nbt.4200
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, yHat] = kalmanFilter(y, A, C, K)

    % Decodes mood state over time from neural features
    % Inputs:
    % - (1) y (time x data): observation timeseries
    % - (2 and later) model parameters, as follows:
    %           - (4) A: LSSM state equation A matrix
    %           - (5) K: LSSM state equation K matrix (Kalman gain)
    %           - (6) C: LSSM observation equation C matrix

    % Outputs:
    % - (1) Xp (time x 1): Predicted states

    N = size(y, 1); % Total number of data points
    
    % Run Kalman filter 
    nx = size(A, 1); % Dimension of hidden state
    x = nan(N, nx);  % Store estimated hidden states here
    
    Xp = zeros(nx, 1);
    for i = 1:N
        x(i, :) = Xp;          % Predicted state x(i|i-1)
        yi = y(i, :)';         % Observation y(i)
        if ~any(isnan(yi))
            zi = yi - C*Xp;        % Innovation z(i)
            Xp = A * Xp + K * zi;  % Prediction of next step
        else
            Xp = A * Xp;  % Prediction of next step
        end
    end
    
    yHat = (C * x')';
end
