function disp_params(params_marginals, initial_cond, flag)
% Disply on the console of the 1st calibration params
% 
% INPUTS
% params_marginals:     [VECTOR]the calibrated parameters for the marginals
% nu_1:                 [SCALAR]the calibrated nu_1
% nu_2:                 [SCALAR]the calibrated nu_2
% nu_z:                 [SCALAR]the calibrated nu_z
% initial_cond:         [VECTOR]the initial condition used
% flag:                 [BOOLEAN]a boolean to save the results to a text file
%
% OUTPUTS
% None
%
% USES      none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    params_USA = params_marginals(1:3);
    params_EU = params_marginals(4:6);

    kappa_USA = params_USA(1);
    theta_USA = params_USA(2);
    sigma_USA = params_USA(3);

    kappa_EU = params_EU(1);
    theta_EU = params_EU(2);
    sigma_EU = params_EU(3);

    if flag == 1

        fileID = fopen('results.txt', 'w');

        % Writing and displaying the initial condition
        fprintf(fileID,'PARAMETERS obtained by calibration:\n');
        fprintf(fileID, 'X0 used :\n');
        fprintf(fileID, '%f \n', initial_cond');
        fprintf(fileID, '-----------------------\n');

        % Writing and displaying the USA market parameters
        fprintf(fileID, 'Calibrated parameters for the USA market: \n');
        fprintf(fileID, 'kappa_USA: %f \n', kappa_USA);
        fprintf(fileID, 'theta_USA: %f \n', theta_USA);
        fprintf(fileID, 'sigma_USA: %f \n', sigma_USA);
        fprintf(fileID, '-----------------------\n');

        % Writing and displaying the EU market parameters
        fprintf(fileID, 'Calibrated parameters for the EU market: \n');
        fprintf(fileID, 'kappa_EU: %f \n', kappa_EU);
        fprintf(fileID, 'theta_EU: %f \n', theta_EU);
        fprintf(fileID, 'sigma_EU: %f \n', sigma_EU);
        fprintf(fileID, '-----------------------\n');


        fclose(fileID);

    end

    % Display initial condition
    disp('-----------------------')        
    disp('X0 used:')
    fprintf('%f\n', initial_cond)

    fprintf('PARAMETERS obtained by calibration:\n');

    disp('-----------------------')
        
    % Display calibrated parameters for the USA market
    disp('Calibrated parameters for the USA market:')
    fprintf('kappa_USA: %f\n', kappa_USA)
    fprintf('theta_USA: %f\n', theta_USA)
    fprintf('sigma_USA: %f\n', sigma_USA)
    disp('-----------------------')
        
    % Display calibrated parameters for the EU market
    disp('Calibrated parameters for the EU market:')
    fprintf('kappa_EU: %f\n', kappa_EU)
    fprintf('theta_EU: %f\n', theta_EU)
    fprintf('sigma_EU: %f\n', sigma_EU)
    disp('-----------------------')

end % function disp_params