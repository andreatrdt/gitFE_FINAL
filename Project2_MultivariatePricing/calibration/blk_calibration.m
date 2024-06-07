function MSE_total = blk_calibration(sigma, dataset, F0, B0, date_settlement)
% Computation of the total RMSE for the dataset given
% 
% INPUT:
% sigma:             [SCALAR] volatility of the process
% dataset:           [STRUCT] dataset preprocessed
% F0:                [VECTOR] forward values
% B0:                [VECTOR] discount values
% date_settlement:   [DATENUM] initial date
% 
% OUTPUT:
% MSE_total:        error to minimize
%
% USES: blkprice()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Conventions
    conv_ACT365 = 3;

    %% Initial parameters
    MSE = zeros(length(dataset.datesExpiry), 1);

    %% Computation 
    
    for ii = 1:length(dataset.datesExpiry)

        %% Initialization
        put_length = length(dataset.putAsk(ii).prices);

        % Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        
        % Interest rate
        interest_rate = -log(B0(ii))/TTM;

        %% Pricing 
        
        % Price of Call/Puts through the Black formula
        [call_prices, put_prices] = blkprice(F0(ii), dataset.strikes(ii).value, interest_rate, TTM, sigma);
        
        % Selection of the OTM options
        call_prices = call_prices(put_length+1:end);
        put_prices = put_prices(1:put_length);
        
        % Mid market prices
        mean_call_price = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        mean_put_price = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;

        %% Computation of RMSE

        MSE(ii) = rmse([call_prices'; put_prices'],[mean_call_price'; mean_put_price']);
    end

    %% Final adjusting of RMSE
    MSE_total = sum(MSE);

end % function blk_calibration