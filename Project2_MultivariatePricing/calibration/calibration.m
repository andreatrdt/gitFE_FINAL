function distance = calibration(params, data_EU, data_USA, F0_EU, B0_EU, F0_USA, B0_USA, date_settlement)
% Computation of the function handle for the calibration procedure
% 
% INPUT:
% params:             [VECTOR] parameters that needs to be calibrated, they're written as
%                     [k1, theta1, sigma1, k2, theta2, sigma2]
% data_EU:            [STRUCT] dataset of the European mkt
% data_USA:           [STRUCT] dataset of the American mkt
% F0_EU:              [VECTOR] initial forward F(0, T) EU
% B0_EU:              [VECTOR] initial discount B(0, T) EU
% F_USA:              [VECTOR] initial forward F(0, T) USA
% B0_USA:             [VECTOR] initial discount B(0, T) USA 
% date_settlement:    [DATENUM] settlement date of the computation
% 
% OUTPUT:
% distance:           [SCALAR] value to be reduced to 0
% 
% USES:              RMSE_total()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Unpacking of the parameters

    params_USA = params(1:3);
    params_EU = params(4:6);

    %% Computation of the RMSE

    RMSE_EU = RMSE_total(params_EU, data_EU, F0_EU, B0_EU, date_settlement);
    RMSE_USA = RMSE_total(params_USA, data_USA, F0_USA, B0_USA, date_settlement); 

    %% Computation of the weights

    weight_USA = data_USA.spot/(data_EU.spot+data_USA.spot);
    weight_EU = data_EU.spot/(data_EU.spot+data_USA.spot);
   
    %% Computation of the distance

    distance = weight_EU*RMSE_EU + weight_USA*RMSE_USA;

end % function calibration