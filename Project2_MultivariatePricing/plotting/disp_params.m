function disp_params(params_marginals, nu_1 ,nu_2 ,nu_z, initial_cond, flag)
%
% INPUTS
% params_marginals: the calibrated parameters for the marginals
% nu_1: the calibrated nu_1
% nu_2: the calibrated nu_2
% nu_z: the calibrated nu_z
% initial_cond: the initial condition used
% flag: a boolean to save the results to a text file
%
% OUTPUTS
% None
%
% USES disp_params()

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
        fprintf(fileID, '-----------------------\n');
        fprintf(fileID, 'X0 used :\n');
        fprintf(fileID, '%f \n', initial_cond);
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

        % Writing and displaying the calibrated nu parameters
        fprintf(fileID, 'calibrated nu_1: \n');
        fprintf(fileID, '%f \n', nu_1);
        fprintf(fileID, '-----------------------\n');
        fprintf(fileID, 'calibrated nu_2: \n');
        fprintf(fileID, '%f \n', nu_2);
        fprintf(fileID, '-----------------------\n');
        fprintf(fileID, 'calibrated nu_z: \n');
        fprintf(fileID, '%f \n', nu_z);
        fprintf(fileID, '-----------------------\n');

        fclose(fileID);

    end

    % Display initial condition
    disp('-----------------------')        
    disp('X0 used:')
    fprintf('%f\n', initial_cond)
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
        
    % Display calibrated nu parameters
    disp('Calibrated nu parameters:')
    fprintf('nu_1: %f\n', nu_1)
    disp('-----------------------')
    fprintf('nu_2: %f\n', nu_2)
    disp('-----------------------')
    fprintf('nu_z: %f\n', nu_z)
    disp('-----------------------')
    


end