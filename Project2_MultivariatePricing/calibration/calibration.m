function calibration(data_EU, data_USA, F_0_EU, discount_at_expiry_EU, date_settlement, F_0_USA, discount_at_expiry_USA , dates_EU, dates_USA, conv_ACT365)

    % EU:
    data = data_EU;

    RMSE_EU = @(p) 0.*p;

    RMSE_USA = @(p) 0.*p;

    for idx = 1:length(dates_EU)

        put_length = length(data.putAsk);
        % compute the log moneyess from the strikes
        log_moneyness = log(F_0_EU(idx) ./ data.strikes(idx).value);

        % time to maturity
        t = yearfrac(date_settlement,dates_EU(idx),conv_ACT365);

        % create a function that the prices of the call options given the strikes

        if idx <= put_length
            prices_EU = @(p) callIntegral(discount_at_expiry_EU(idx), F_0_EU(idx), p(1), p(2), p(3), t, log_moneyness)-data.strikes(idx).value + data.spot; 
        else
            prices_EU = @(p) callIntegral(discount_at_expiry_EU(idx), F_0_EU(idx), p(1), p(2), p(3), t, log_moneyness);
        end

        mean_call_price_EU = (data.callAsk(idx).prices+data.callBid(idx).prices)/2;

        RMSE_EU = @(p) RMSE_EU(p) + sum((prices_EU(p) - mean_call_price_EU).^2);
    end


    RMSE_EU =  @(p) sqrt(RMSE_EU(p))

    % USA:

    data = data_USA;

    for idx = 1:length(dates_USA)

        put_length = length(data.putAsk);

        % compute the log moneyess from the strikes
        log_moneyness = log(F_0_USA(idx) ./ data.strikes(idx).value);

        % time to maturity
        t = yearfrac(date_settlement,dates_USA(idx),conv_ACT365);

        % create a function that the prices of the call options given the strikes
        if idx <= put_length
            prices_USA = @(p) callIntegral(discount_at_expiry_USA(idx), F_0_USA(idx), p(1), p(2), p(3), t, log_moneyness)-data.strikes(idx).value + data.spot; 
        else
            prices_USA = @(p) callIntegral(discount_at_expiry_USA(idx), F_0_USA(idx), p(1), p(2), p(3), t, log_moneyness);
        end
        
        mean_call_price_USA = (data_USA.callAsk(idx).prices+data_USA.callBid(idx).prices)/2;

        RMSE_USA = @(p) RMSE_USA(p) + sum((prices_USA(p) - mean_call_price_USA).^2);
    end

    RMSE_USA =  @(p) sqrt(RMSE_USA(p))


    weights = [data_EU.spot/(data_EU.spot+data_USA.spot), data_USA.spot/(data_EU.spot+data_USA.spot)];

    dist(p) = weights(1)*RMSE_EU(p) + weights(2)*RMSE_USA(p);

    % Calibration
    x0 = ones(6, 1);
 
    % Linear inequality constraints on the theta_i
    A = []; b = [];
 
    % Unused inequality matrixies
    Aeq = []; beq = [];

    % Bounds for the single parameter, no ub required
    lb = [0; -Inf; 0; 0; -Inf; 0];
    ub = [];

    % Options for the visualization
    options = optimset('Display', 'iter');

    x = fmincon(@(x) dist(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconstr(x), options);

end