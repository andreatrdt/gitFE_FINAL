function dataset_exploration(data_EU, data_USA)
% Exploration of the dataset and consequent plots
% 
% INPUT:
% data_EU:               [STRUCT]dataset EU
% data_USA:              [STRUCT]dataset USA
%
% OUTPUT:
% None
%
% USES:         none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Plot of the European Call Prices

    figure();

    for ii = 1:length(data_EU.datesExpiry)

        % Compute the mid prices
        mid_call = (data_EU.callAsk(ii).prices + data_EU.callBid(ii).prices)/2;
        plot(data_EU.strikes(ii).value, mid_call); hold on;
    end

    grid on; title('European Call Prices');
    xlabel('Strikes'); ylabel('Prices');

    %% Plot of the European Put Prices

    figure();

    for ii = 1:length(data_EU.datesExpiry)

        % Compute the mid prices
        mid_put = (data_EU.putAsk(ii).prices + data_EU.putBid(ii).prices)/2;
        plot(data_EU.strikes(ii).value, mid_put); hold on;
    end

    grid on; title('European Put Prices');
    xlabel('Strikes'); ylabel('Prices');

    %% Plot of the American Call Prices

    figure();

    for ii = 1:length(data_USA.datesExpiry)

        % Compute the mid prices
        mid_call = (data_USA.callAsk(ii).prices + data_USA.callBid(ii).prices)/2;
        plot(data_USA.strikes(ii).value, mid_call); hold on;
    end

    grid on; title('American Call Prices');
    xlabel('Strikes'); ylabel('Prices');

    %% Plot of the American Put Prices

    figure();

    for ii = 1:length(data_USA.datesExpiry)

        % Compute the mid prices
        mid_put = (data_USA.putAsk(ii).prices + data_USA.putBid(ii).prices)/2;
        plot(data_USA.strikes(ii).value, mid_put); hold on;
    end

    grid on; title('American Put Prices');
    xlabel('Strikes'); ylabel('Prices');

    %%

    %% Study of the penny options
    % Penny options are those options with value lower than 0.1 index
    % points
    
    % European

    count_EU_call_ask = 0;
    count_EU_call_bid = 0;
    count_EU_put_ask = 0;
    count_EU_put_bid = 0;

    for ii = 1:length(data_EU.datesExpiry)

        count_EU_call_ask = count_EU_call_ask + sum(data_EU.callAsk(ii).prices <= 0.1);

        count_EU_call_bid = count_EU_call_bid + sum(data_EU.callBid(ii).prices <= 0.1);

        count_EU_put_ask = count_EU_put_ask + sum(data_EU.putAsk(ii).prices <= 0.1);

        count_EU_put_bid = count_EU_put_bid + sum(data_EU.putBid(ii).prices <= 0.1);
    end

    disp('The number of European Penny options is');
    disp(count_EU_call_ask+count_EU_call_bid+count_EU_put_ask+count_EU_put_bid);

    % American

    count_USA_call_ask = 0;
    count_USA_call_bid = 0;
    count_USA_put_ask = 0;
    count_USA_put_bid = 0;

    for ii = 1:length(data_USA.datesExpiry)

        count_USA_call_ask = count_USA_call_ask + sum(data_USA.callAsk(ii).prices <= 0.1);

        count_USA_call_bid = count_USA_call_bid + sum(data_USA.callBid(ii).prices <= 0.1);

        count_USA_put_ask = count_USA_put_ask + sum(data_USA.putAsk(ii).prices <= 0.1);

        count_USA_put_bid = count_USA_put_bid + sum(data_USA.putBid(ii).prices <= 0.1);
    end

    disp('The number of American Penny options is');
    disp(count_USA_call_ask+count_USA_call_bid+count_USA_put_ask+count_USA_put_bid);

    %% Study of the liquidity criterion
    % The liquidity criterion ask for the bid/ask difference to be lower
    % than the ratio (ask - bid)/ask

    % European

    indicator_liq_call = 0;
    indicator_liq_put = 0;

    for ii = 1:length(data_EU.datesExpiry)

        liquidity_call = (data_EU.callAsk(ii).prices - data_EU.callBid(ii).prices)/data_EU.callAsk(ii).prices;
        indicator_liq_call = indicator_liq_call + sum(liquidity_call >= 0.6);

        liquidity_put = (data_EU.putAsk(ii).prices - data_EU.putBid(ii).prices)/data_EU.putAsk(ii).prices;
        indicator_liq_put = indicator_liq_put + sum(liquidity_put >= 0.6);

    end

    disp('The number of illiquid options in the EU mkt is:');
    disp(indicator_liq_call + indicator_liq_put);

    % American
    
    indicator_liq_call = 0;
    indicator_liq_put = 0;

    for ii = 1:length(data_USA.datesExpiry)

        liquidity_call = (data_USA.callAsk(ii).prices - data_USA.callBid(ii).prices)/data_USA.callAsk(ii).prices;
        indicator_liq_call = indicator_liq_call + sum(liquidity_call >= 0.6);

        liquidity_put = (data_USA.putAsk(ii).prices - data_USA.putBid(ii).prices)/data_USA.putAsk(ii).prices;
        indicator_liq_put = indicator_liq_put + sum(liquidity_put >= 0.6);

    end

    disp('The number of illiquid options in the USA mkt is:');
    disp(indicator_liq_call + indicator_liq_put);
    
end % function dataset_exploration