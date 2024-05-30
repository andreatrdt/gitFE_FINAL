function dataset_exploration(data_EU, data_USA, date_settlement)
% Exploration of the dataset and consequent plots
% 
% INPUT:
% data_EU:               dataset EU
% data_USA:              dataset USA
% date_settlement:       initial date of the computation

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

end