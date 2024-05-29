function disp_params(params_marginals, nu_1 ,nu_2 ,nu_z, initial_cond, flag)

    params_USA = params_marginals(1:3);
    params_EU = params_marginals(4:6);

    kappa_USA = params_USA(1);
    theta_USA = params_USA(2);
    sigma_USA = params_USA(3);


    kappa_EU = params_EU(1);
    theta_EU = params_EU(2);
    sigma_EU = params_EU(3);

    if flag == 1

        fileID = fopen('results.txt', 'a');

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

    fprintf('-----------------------\n');
    fprintf('X0 used :\n');
    fprintf('%f \n', initial_cond);
    fprintf('-----------------------\n');

    fprintf('Calibrated parameters for the USA market: \n');
    fprintf('kappa_USA: %f \n', kappa_USA);
    fprintf('theta_USA: %f \n', theta_USA);
    fprintf('sigma_USA: %f \n', sigma_USA);
    fprintf('-----------------------\n');

    
    fprintf('Calibrated parameters for the EU market: \n');
    fprintf('kappa_EU: %f \n', kappa_EU);
    fprintf('theta_EU: %f \n', theta_EU);
    fprintf('sigma_EU: %f \n', sigma_EU);
    fprintf('-----------------------\n');

    fprintf('calibrated nu_1: \n');
    fprintf('%f \n', nu_1);
    fprintf('-----------------------\n');
    fprintf('calibrated nu_2: \n');
    fprintf('%f \n', nu_2);
    fprintf('-----------------------\n');
    fprintf('calibrated nu_z: \n');
    fprintf('%f \n', nu_z);
    fprintf('-----------------------\n');


end