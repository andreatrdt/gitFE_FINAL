function disp_contract_prices(price_Levy,CI_Levy,price_Levy_AV, CI_Levy_AV, ...
    price_Blk,CI_Blk,price_Blk_AV,CI_Blk_AV, price_Levy_closed, price_SemiclosedBlk)
%
% INPUT:
% price_Levy:           [SCALAR] price of the derivative under the Levy model
% CI_Levy:              [VECTOR] confidence interval for the Levy model
% price_Blk:            [SCALAR] price of the derivative under the Black model
% CI_Blk:               [VECTOR] confidence interval for the Black model
% price_Blk_AV:         [SCALAR] price of the derivative under the Black model with AV
% CI_Blk_AV:            [VECTOR] confidence interval for the Black model with AV
% price_SemiclosedBlk:  [SCALAR] price of the derivative under the Black model with semi-closed formula
%
% OUTPUT:
% Display the prices of the derivative and their confidence intervals
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    disp('Prices of the derivative')
    disp('-----------------------------------------------------------------------------------------------------')
    fprintf(' Model              |  Price              |  Confidence Interval \n')
    disp('-----------------------------------------------------------------------------------------------------')
    fprintf('Levy:               |  %.8f        |  [ %.8f ,  %.8f] \n', price_Levy, CI_Levy(1), CI_Levy(2))

    fprintf('Levy AV:            |  %.8f        |  [ %.8f ,  %.8f] \n', price_Levy_AV, CI_Levy_AV(1), CI_Levy_AV(2))
    
    fprintf('Levy semi-closed:   |  %.8f        |  -- \n', price_Levy_closed)

    fprintf('Black:              |  %.8f        |  [ %.8f ,  %.8f] \n', price_Blk, CI_Blk(1), CI_Blk(2))
    
    fprintf('Black AV:           |  %.8f        |  [ %.8f ,  %.8f] \n', price_Blk_AV, CI_Blk_AV(1), CI_Blk_AV(2))
    
    fprintf('Black semi-closed:  |  %.8f        |  -- \n', price_SemiclosedBlk)

    disp('-----------------------------------------------------------------------------------------------------')    

end % function disp_contract_prices