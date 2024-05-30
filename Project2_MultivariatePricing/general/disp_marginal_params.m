function disp_marginal_params(sol_USA , sol_EU ,Beta_z, gamma_z, save_results)
%
% INPUTS
% sol_USA: structure containing the solution of the optimization problem
% sol_EU: structure containing the solution of the optimization problem
% nu_z: the calibrated nu_z
% save_results: a boolean to save the results to a text file
%
% OUTPUTS
% None
%
% USES disp_marginal_params()


     % Extract parameters for the USA
     a_USA = sol_USA(1);
     Beta_USA = sol_USA(2);
     gamma_USA = sol_USA(3);
  
     % Extract parameters for the EU
     a_EU = sol_EU(1);
     Beta_EU = sol_EU(2);
     gamma_EU = sol_EU(3);



    % Display the parameters

    fprintf('----------------------------\n');
    fprintf('a_USA: %.4f\n', a_USA);
    fprintf('a_EU: %.4f\n', a_EU);
    fprintf('----------------------------\n');
    fprintf('Beta_z: %.4f\n', Beta_z);
    fprintf('Beta_USA: %.4f\n', Beta_USA);
    fprintf('Beta_EU: %.4f\n', Beta_EU);
    fprintf('----------------------------\n');
    fprintf('gamma_z: %.4f\n', gamma_z);
    fprintf('gamma_USA: %.4f\n', gamma_USA);
    fprintf('gamma_EU: %.4f\n', gamma_EU);
    fprintf('----------------------------\n');

    if save_results
        % Save the parameters to a text file
        fileID = fopen('results.txt', 'a');

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

        % Close the file
        fclose(fileID);
    end
end