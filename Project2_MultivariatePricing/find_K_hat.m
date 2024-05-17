function K_hat = find_K_hat(strikes,index)
% 
% This function finds K_hat
%
% Input:
% data_EU: data for the European options
% data_USA: data for the American options
% Output:
% K_hat

% Find the synthetic forward hat operator G_hat

    K_hat = sum(strikes(index).value)/length(strikes(index).value);

end % function find_K_hat