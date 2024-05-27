function rho = hist_corr(dataset)
% Computation of the historical correlation
% 
% INPUT:
% dataset:        [STRUCT] initial dataset of the returns
% 
% OUTPUT:
% rho:            correlation coefficient between USA and EU

    %% Computation of the returns

    returns_USA = dataset.Annually(:,1);
    returns_EU = dataset.Annually(:,2);
    
    %% Computation of the correlation coefficient

    rho = corr(returns_EU,returns_USA);

end % function hist_corr