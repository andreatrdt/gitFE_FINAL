function plot_returns(data,date_settlement)
% Plotting of the data returns from the S&P500 and EUROSTOXX50 dataset
% 
% INPUTS:
% data:                 [STRUCT] struct od the DATA used
% date_settlement:      [DATENUM] the date of the settlement
% 
% OUTPUTS: none
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Load the data
    returns_EU = flip(data.Daily(:,1));

    returns_USA = flip(data.Daily(:,2));

    % Create the dates
    days = [date_settlement : date_settlement+length(returns_USA)-1];

    %% Plot the returns
    figure
    subplot(3,1,1)
    plot(days,returns_EU)
    title('EU returns')
    datetick('x','dd-mmm-yyyy','keepticks')

    subplot(3,1,2)
    plot(days,returns_USA,"Color","r")
    title('USA returns')
    datetick('x','dd-mmm-yyyy','keepticks')
    
    subplot(3,1,3)
    plot(days,returns_EU); hold on;
    plot(days,returns_USA)
    title('EU and USA returns')
    legend('EU returns', 'USA returns')
    datetick('x','dd-mmm-yyyy','keepticks')
    hold off
    

end % function plot_returns