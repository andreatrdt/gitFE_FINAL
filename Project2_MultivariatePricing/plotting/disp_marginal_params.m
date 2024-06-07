function disp_marginal_params(sol_USA , sol_EU ,Beta_z, gamma_z,nu_z, save_results)
%
% INPUTS
% sol_USA:          [STRUCT]structure containing the solution of the optimization problem
% sol_EU:           [STRUCT]structure containing the solution of the optimization problem
% Beta_z:           [SCALAR]the calibrated beta_z
% gamma_z:          [SCALAR]the calibrated gamma
% save_results:     [BOOLEAN]a boolean to save the results to a text file
%
% OUTPUTS
% None
%
% USES none

% Authors:
% M.Maspes, A.Tarditi, M.Torba




    % Extract parameters for the USA
    a_USA = sol_USA(4);
    Beta_USA = sol_USA(2);
    gamma_USA = sol_USA(3);
    nu_1 = sol_USA(1);

    % Extract parameters for the EU
    a_EU = sol_EU(4);
    Beta_EU = sol_EU(2);
    gamma_EU = sol_EU(3);
    nu_2 = sol_EU(1);

    % Display the parameters
    fprintf('PARAMETERS obtained by convolution:\n')
    disp('----------------------------')
    fprintf('a_USA: %.4f\n', a_USA)
    fprintf('a_EU: %.4f\n', a_EU)
    disp('----------------------------')
    
    % Display Beta parameters
    fprintf('Beta_z: %.4f\n', Beta_z)
    fprintf('Beta_USA: %.4f\n', Beta_USA)
    fprintf('Beta_EU: %.4f\n', Beta_EU)
    disp('----------------------------')
    
    % Display gamma parameters
    fprintf('gamma_z: %.4f\n', gamma_z)
    fprintf('gamma_USA: %.4f\n', gamma_USA)
    fprintf('gamma_EU: %.4f\n', gamma_EU)
    disp('----------------------------')

    fprintf('nu_1: %.4f \n', nu_1);
    fprintf('nu_2: %.4f \n', nu_2);
    fprintf('nu_z: %.4f \n', nu_z);
    disp('----------------------------')


    if save_results
        % Save the parameters to a text file
        fileID = fopen('results.txt', 'a');
        
        fprintf(fileID,'PARAMETERS obtained by convolution:\n');
        fprintf(fileID,'----------------------------\n');
        fprintf(fileID,'a_USA: %.4f\n', a_USA);
        fprintf(fileID,'a_EU: %.4f\n', a_EU);
        fprintf(fileID,'----------------------------\n');
        fprintf(fileID,'Beta_z: %.4f\n', Beta_z);
        fprintf(fileID,'Beta_USA: %.4f\n', Beta_USA);
        fprintf(fileID,'Beta_EU: %.4f\n', Beta_EU);
        fprintf(fileID,'----------------------------\n');
        fprintf(fileID,'gamma_z: %.4f\n', gamma_z);
        fprintf(fileID,'gamma_USA: %.4f\n', gamma_USA);
        fprintf(fileID,'gamma_EU: %.4f\n', gamma_EU);
        fprintf(fileID,'----------------------------\n');
        fprintf(fileID, 'nu_1: %.4f \n', nu_1);
        fprintf(fileID, 'nu_2: %.4f \n', nu_2);
        fprintf(fileID, 'nu_z: %.4f \n', nu_z);
        fprintf(fileID,'----------------------------\n');

        % Close the file
        fclose(fileID);
    end

end % function disp_marginal_params