function plot_multi_calib(dates_EU, dates_USA, fun_eval, length_EU, length_USA)
% Plot of the evaluation errors of the multicalibration procedure
% 
% INPUT:
% dates_EU:        dates expiry of EU
% dates_USA:       dates expiry of USA
% fun_eval:        eval obj func from multicalib
% length_EU:       # european maturities
% length_USA:      # american maturities
% 
% OUTPUT: none
% 
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    figure;
    plot(dates_EU, fun_eval(1 : length_EU), '*-'); grid on;
    datetick('x','dd-mmm-yyyy','keepticks');
    title('EU - Obj func minimum');
    ylabel('Obj function value');

    figure;
    plot(dates_USA(1:length_USA), fun_eval(length_EU+1 : length_USA+length_EU), '*-'); grid on;
    datetick('x','dd-mmm-yyyy','keepticks');
    title('USA Obj func minimum');
    ylabel('Obj function value');

end % function plot_multi_calib