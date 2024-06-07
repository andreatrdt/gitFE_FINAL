function surface_vols(data,F0)
% Plot surface volatilities in 3D for better comprehension of the
% distributions
%
% INPUTS:
% data    [STRUCT]  structure containing the data
% F0      [SCALAR]  forward price
%
% OUTPUTS: none
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    %% Initialize vector to store the number of strikes for each expiry date
    strikes = [];

    % Consolidate all unique strikes across all expiry dates
    for i = 1:length(data.datesExpiry)
        strikes = unique([strikes, data.strikes(i).value]);
    end

    % Sort the strikes array
    strikes = sort(strikes);

    %% Initialize the complete volatility matrix with NaNs
    vols = NaN(length(data.datesExpiry), length(strikes));

    index_matrix = NaN(length(data.datesExpiry), length(strikes));

    % Fill the volatility matrix with the corresponding volatilities
    for i = 1:length(data.datesExpiry)
        % Find the indices of current strikes in the consolidated strikes array
        [~, idxes] = ismember(data.strikes(i).value, strikes);

        % Combine putBid and callAsk implied volatilities
        volatility = [data.putBid(i).impvol, data.callAsk(i).impvol];

        % Populate the vols matrix with the corresponding volatilities
        index_matrix(i, idxes) = idxes; 
        vols(i, idxes) = volatility;
    end


    %% Perform 2D interpolation to fill missing values
    [X, Y] = meshgrid(strikes, datenum(data.datesExpiry));
    validMask = ~isnan(vols);
    interpolatedVols = griddata(X(validMask), Y(validMask), vols(validMask), X, Y, 'linear');

    log_moneyness = log(F0./strikes);

    % Convert dates to strings for the y-axis labels
    dateStrings = data.datesExpiry;

    %% Plot the surface
    figure();
    surf(log_moneyness, datenum(data.datesExpiry),interpolatedVols,'FaceAlpha',0.5);
    hold on
    mesh(log_moneyness, datenum(data.datesExpiry),interpolatedVols);
    hold off
    
    % Label the axes
    xlabel('Strikes');
    ylabel('Dates');
    zlabel('Implied Volatility');
    title('Surface of Implied Volatilities');

    % Set the y-axis ticks and labels
    datetick('y','dd-mmm-yyyy','keepticks')
    colormap("parula")

end % function surface_vols()

