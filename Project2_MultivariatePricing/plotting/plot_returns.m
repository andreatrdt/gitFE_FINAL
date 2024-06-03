function plot_returns(data,date_settlement)
% Plotting of the data returns
% 
% INPUTS:
% data: struct od the DATA used
% date_settlement: the date of the settlement
% 
% OUTPUTS:
% plot the returns of the options

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

    % filename_EU = 'eurostocks.xlsx';

    % filename_USA = 'SP500.xls';

    % returns_EU = readtable(filename_EU,Range="C6:C229");

    % returns_EU = table2array(returns_EU);

    % dates_EU = readtable(filename_EU,RAnge = "B6:B229");

    % dates_EU = table2array(dates_EU);

  
    % returns_USA = readtable(filename_USA,Range = "B12:B1317");

    % returns_USA = table2array(returns_USA);

    % dates_USA = readtable(filename_USA,Range = "A12:A1317");

    % dates_USA = table2array(dates_USA);

    
    % figure
    % subplot(2,1,1)
    % plot(datenum(dates_EU),returns_EU)
    % title('EU real returns')
    % datetick('x','dd-mmm-yyyy','keepticks');
    % subplot(2,1,2)
    % plot(datenum(dates_USA),returns_USA,"Color","r")
    % title('USA real returns')
    % datetick('x','dd-mmm-yyyy','keepticks');
    


end % function plot_returns