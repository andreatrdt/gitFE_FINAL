function rho = hist_corr(dataset)
% Computation of the historical correlation through the Annualy returns on
% the S&P500 and EUROSTOXX50 dataset
% 
% INPUT:
% dataset:        [STRUCT] initial dataset of the returns
% 
% OUTPUT:
% rho:            [SCALAR] correlation coefficient between USA and EU
%
% USES:        none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Computation of the returns

    returns_USA = dataset.Annually(:,1);
    returns_EU = dataset.Annually(:,2);
    
    %% Computation of the correlation coefficient

    rho = corr(returns_EU,returns_USA);

end % function hist_corr