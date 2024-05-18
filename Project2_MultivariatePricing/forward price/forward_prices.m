function [F_vector, G_vector] = forward_prices(dataset, date)
% Computation of the forward prices following the Baviera, Azzone paper
% 
% INPUT:
% dataset:           data containing all the required tables
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

    %% Computation of the strikes and synthetic forwards

    % Strikes
    Ki = dataset.strikes(index).value;

    % Synthetic forwards
    [Gi, Gi_ask, Gi_bid] = synthethic_forward(dataset.callBid, dataset.callAsk, ...
        dataset.putBid, dataset.putAsk, index);
    
    %% Computation of the estimated discount factor

    B_bar = estimation_discount_factor(Gi, Ki);
    
    %% Computation of the forward prices
    F = Gi/B_bar + Ki;
    F_ask = Gi_ask/B_bar + Ki;
    F_bid = Gi_bid/B_bar + Ki;
    
    %% Plot of the figures

%     figure();
%     plot(Ki, F_ask, '*-'); hold on;
%     plot(Ki, F, 'o-');
%     plot(Ki, F_bid, '*-'); 
% 
%     % title('Forward prices at Expiry ', date);
%     xlabel('Strikes'); ylabel('Forwards prices');
%     legend('Ask', 'Mid', 'Bid')
%     grid on; hold off;

    %% Computation of the slope

    p = polyfit(Ki, F - Ki, 1);
    slope = p(1);
    disp(['Estimated Slope: ', num2str(slope)]);

    %% Creation vectors to return
    F_vector = [F; F_ask; F_bid];
    G_vector = [Gi; Gi_ask; Gi_bid];

end % function forward_prices