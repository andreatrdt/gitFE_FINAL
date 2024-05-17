function G_hat = find_synthetic_forward_hat(G)
% 
% This function finds the synthetic forward hat operator G_hat
%
% INPUT:
% data_EU: data for the European options
% data_USA: data for the American options
% OUTPUT:
% G_hat: synthetic forward hat operator

% Find the synthetic forward hat operator G_hat

    G_hat = sum(G)/length(G);

end % function find_synthetic_forward_hat