function [F_vector, B_bar_vector] = forward_prices(dataset, flag)
% Computation of the forward prices following the Baviera, Azzone paper
% 
% INPUT:
% dataset:           [STRUCT] data containing all the required tables
% date:              [SCALAR] date to consider for the computation
% flag:              [0: without plots & slope; 1: with plots & slope]
% 
% OUTPUT:
% F_vector:          [VECTOR] forwards value F(0, T)
% B_bar_vector:      [VECTOR] market calibrated discount B(0, T)
% 
% USES:     synthethic_forward(), estimation_discount_factor()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Introduction of the return vectors

    F_vector = zeros(length(dataset.datesExpiry), 1);
    B_bar_vector = zeros(length(dataset.datesExpiry), 1);

    %% Computation of the forwards

    for ii = 1:length(dataset.datesExpiry)

        %% Computation of the strikes and synthetic forwards

        % Strikes
        Ki = dataset.strikes(ii).value;
    
        % Synthetic forwards
        [Gi, Gi_ask, Gi_bid] = synthethic_forward(dataset.callBid(ii).prices, dataset.callAsk(ii).prices, ...
            dataset.putBid(ii).prices, dataset.putAsk(ii).prices);

        %% Computation of the estimated discount factor

        B_bar_vector(ii) = estimation_discount_factor(Gi, Ki);

        %% Computation of the forward prices
        F_i_vector = Gi/B_bar_vector(ii) + Ki;
        F_ask__i_vector = Gi_ask/B_bar_vector(ii) + Ki;
        F_bid__i_vector = Gi_bid/B_bar_vector(ii) + Ki;
        
        % Computation of the required Forwards
        F_vector(ii) = mean(F_i_vector);

        %% Eventual plot

        if flag
            figure();
            plot(Ki, F_ask__i_vector, '*-'); hold on;
            plot(Ki, F_i_vector, 'o-');
            plot(Ki, F_bid__i_vector, '*-'); grid on;
        
            xlabel('Strikes'); ylabel('Forwards prices');
            legend('Ask', 'Mid', 'Bid')
            grid on; hold off;
        end
        
        %% Computation of the slope

        if flag
            p = polyfit(Ki, F_i_vector - Ki, 1);
            slope = p(1);
            disp(['Estimated Slope: ', num2str(slope)]);
        end

    end

end % function forward_prices