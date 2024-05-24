function surface_vols(data)
    % Initialize vector to store the number of strikes for each expiry date
    vect = zeros(length(data.datesExpiry), 1);

    % Find the maximum number of strikes for any expiry date
    for i = 1:length(data.datesExpiry)
        vect(i) = length(data.strikes(i).value);
    end

    % Identify the expiry date with the maximum number of strikes
    [~, idx] = max(vect);

    % Use the strikes from the expiry date with the maximum number of strikes
    strikes = data.strikes(idx).value;

    % Initialize the volatility matrix with NaNs
    vols = NaN(length(data.datesExpiry), length(strikes));

    % Consolidate all unique strikes across all expiry dates
    for i = 1:length(data.datesExpiry)
        strikes = unique([strikes, data.strikes(i).value]);
    end

    % Sort the strikes array
    strikes = sort(strikes);

    % Initialize the complete volatility matrix with NaNs
    vols = NaN(length(data.datesExpiry), length(strikes));

    % Fill the volatility matrix with the corresponding volatilities
    for i = 1:length(data.datesExpiry)
        % Find the indices of current strikes in the consolidated strikes array
        [~, idxes] = ismember(data.strikes(i).value, strikes);

        % Combine putBid and callAsk implied volatilities
        volatility = [data.putBid(i).impvol, data.callAsk(i).impvol];

        % Populate the vols matrix with the corresponding volatilities
        vols(i, idxes) = volatility;
    end


    % Perform 2D interpolation to fill missing values
    [X, Y] = meshgrid(strikes, datenum(data.datesExpiry));
    validMask = ~isnan(vols);
    interpolatedVols = griddata(X(validMask), Y(validMask), vols(validMask), X, Y, 'linear');



    % Plot the surface
    figure();
    surf(strikes, datenum(data.datesExpiry), interpolatedVols);

    % Label the axes
    xlabel('Strikes');
    ylabel('Dates');
    zlabel('Implied Volatility');
    title('Surface of Implied Volatilities');

    grid on;

    
    shading interp;
    colorbar;
end
