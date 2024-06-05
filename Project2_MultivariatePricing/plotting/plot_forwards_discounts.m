function plot_forwards_discounts(data_EU, data_USA, B_EU, B_USA, F0_EU, F0_USA, date_settlement)
% Showing of the discounts and forwards fot the Latex
% 
% INPUT:
% data_EU:            [STRUCT] data of EU mkt
% data_USA:           [STRUCT] data of USA mkt
% B_EU:               [VECTOR] discounts EU
% B_USA:              [VECTOR] discounts USA
% F0_EU:              [VECTOR] initial forwards EU
% F0_USA.             [VECTOR] initial forwards USA
% date_settlement:    [DATENUM] initial date
%
% OUTPUT:
% None
%
% USES:     none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Conventions

    conv_ACT365 = 3;

    %% Computation of the dates

    dates_EU = datenum(data_EU.datesExpiry);
    dates_USA = datenum(data_USA.datesExpiry);

    %% Plot of the discounts

    figure;
    plot(dates_EU,B_EU,'-*','Color','b'); grid on;
    title('Discount factor for the EU market');
    ylabel('Discounts');
    datetick('x','dd-mmm-yyyy','keepticks')


    figure;
    plot(dates_USA,B_USA,'-*','Color','r'); grid on;
    title('Discount factor for the USA market');
    ylabel('Discounts');
    datetick('x','dd-mmm-yyyy','keepticks')

    %% Plot of the forwards from Synthetic Forwards

    figure;
    plot(dates_EU, F0_EU,'-*', 'Color', 'b'); grid on;
    title('Forwards value EU');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    figure;
    plot(dates_USA, F0_USA, '-*', 'Color', 'b'); grid on;
    title('Forwards value USA');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    %% Plot of the forwards with comparison of the theoretical ones

    figure;
    plot(dates_EU, F0_EU,'-*', 'Color', 'b'); grid on;
    title('Forwards value EU');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    yf = yearfrac(date_settlement, data_EU.datesExpiry, conv_ACT365);
    rates = -log(B_EU)./yf;
    F0_sim_EU = data_EU.spot .* exp(rates .* yf);
    hold on;  plot(dates_EU, F0_sim_EU, '-o', 'Color', 'r');

    figure;
    plot(dates_USA, F0_USA, '-*', 'Color', 'b'); grid on;
    title('Forwards value USA');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    yf = yearfrac(date_settlement, data_USA.datesExpiry, conv_ACT365);
    rates = -log(B_USA)./yf;
    F0_sim_USA = data_USA.spot .* exp(rates .* yf);
    hold on;  plot(dates_USA, F0_sim_USA, '-o', 'Color', 'r');

end % function plot_forwards_discounts