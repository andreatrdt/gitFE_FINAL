function [F_vector, G_vector , B_bar] = forward_prices(dataset, date, flag)
% Computation of the forward prices following the Baviera, Azzone paper
% 
% INPUT:
% dataset:           data containing all the required tables
% date:              date to consider for the computation
% flag:              [0: with plots & slope; 1: without plots & slope]
% 
% OUTPUT:
% F_vector:          [MATRIX] data regarding forward prices
% G_vector:          [MATRIX] data regarding synthetic forward prices
% 
% USES:
% function synthethic_forward()
% function 

    %% Computation of the date index
    % Passing a date from input, we ghet the related value of the index

    index = find_idx(dataset.datesExpiry, date);

    %% Control about penny options and liquidity

    % Necessary to create

    %% Computation of the strikes and synthetic forwards

    % Strikes
    Ki = dataset.strikes(index).value;

    % Synthetic forwards
    [Gi, Gi_ask, Gi_bid] = synthethic_forward(dataset.callBid, dataset.callAsk, ...
        dataset.putBid, dataset.putAsk, index);
    
    %% Computation of the estimated discount factor

    B_bar = estimation_discount_factor(Gi, Ki);
    
    %% Computation of the forward prices
    F_vector = Gi/B_bar + Ki;
    F_ask_vector = Gi_ask/B_bar + Ki;
    F_bid_vector = Gi_bid/B_bar + Ki;
    
    % We're required to compute the minimum and the maximum in order to
    % find the only F0 for the maturity
    F_max = min(F_ask_vector);
    F_min = max(F_bid_vector);

    F = (F_min + F_max)/2;

    %% Plot of the figures
    
    if ~flag
        figure();
        plot(Ki, F_ask_vector, '*-'); hold on;
        plot(Ki, F_vector, 'o-');
        plot(Ki, F_bid_vector, '*-'); grid on;
    
        % title('Forward prices at Expiry ', date);
        xlabel('Strikes'); ylabel('Forwards prices');
        legend('Ask', 'Mid', 'Bid')
        grid on; hold off;
    end

    %% Computation of the slope

    if ~flag
        p = polyfit(Ki, F_vector - Ki, 1);
        slope = p(1);
        disp(['Estimated Slope: ', num2str(slope)]);
    end

    %% Creation vectors to return
    F_vector = [F_vector; F_ask_vector; F_bid_vector];
    G_vector = [Gi; Gi_ask; Gi_bid];

end % function forward_prices