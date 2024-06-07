function B_bar = estimation_discount_factor(Gi, Ki)
% Computation of the estimated discount factor B_bar(0, T)
% 
% INPUT:
% Gi:               [VECTOR] synthetic forwards for the estimation
% Ki:               [VECTOR] strikes for the estimation
% 
% OUTPUT:
% B_bar:            [SCALAR] estimated discount factor for the forward computation
%
% USES:     none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Computation of the mean values
    % Mean values of the strikes and synthetic forwards

    G_hat = mean(Gi);
    K_hat = mean(Ki);

    %% Final computation of B_bar(0, T)
    B_bar =  - (Ki - K_hat) * (Gi - G_hat)' / ((Ki - K_hat) * (Ki - K_hat)');

end % function estimation_discount_factor