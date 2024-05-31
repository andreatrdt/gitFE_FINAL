function disp_contract_prices(price_Levy,CI_Levy,price_Blk,CI_Blk,price_Blk_AV,CI_Blk_AV,price_SemiclosedBlk)
    
    disp('Prices of the derivative')
    disp('-----------------------------------------------------------------------------------------------------')
    fprintf(' Model              |  Price               |  Confidence Interval \n')
    disp('-----------------------------------------------------------------------------------------------------')
    fprintf('Levy:               |  %.8f%%              |  [ %.8f%% ,  %.8f%%] \n', price_Levy, CI_Levy(1), CI_Levy(2))
    
    fprintf('Black:              |  %.8f%%              |  [ %.8f%% ,  %.8f%%] \n', price_Blk, CI_Blk(1), CI_Blk(2))
    
    fprintf('Black AV:           |  %.8f%%              |  [ %.8f%% ,  %.8f%%] \n', price_Blk_AV, CI_Blk_AV(1), CI_Blk_AV(2))
    
    fprintf('Black semi-closed:  |  %.8f%%               |  -- \n', price_SemiclosedBlk)

    disp('-----------------------------------------------------------------------------------------------------')    

end