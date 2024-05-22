function plot_returns(data,date_settlement)
%
% INPUTS:
% data: struct od the DATA used

% OUTPUTS:
% plot the returns of the options

    returns_EU = flip(data.Daily(:,1));

    returns_USA = flip(data.Daily(:,2));

    days = [date_settlement : date_settlement+length(returns_USA)-1];


    figure
    subplot(3,1,1)
    plot(days,returns_EU)
    title('EU returns')
    datetick('x','dd-mmm-yyyy','keepticks')
    subplot(3,1,2)
    plot(days,returns_USA)
    title('USA returns')
    datetick('x','dd-mmm-yyyy','keepticks')
    subplot(3,1,3)
    plot(days,returns_EU)
    hold on
    plot(days,returns_USA)
    title('EU and USA returns')
    legend('EU returns', 'USA returns')
    datetick('x','dd-mmm-yyyy','keepticks')
    hold off




end % function plot_returns