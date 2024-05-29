function disp_marginal_params(sol_USA, sol_EU)
   


    % Extract parameters for the USA
    a_USA = sol_USA.x(1);
    Beta_USA = sol_USA.x(4);
    gamma_USA = sol_USA.x(5);


    % Extract parameters for the EU
    a_EU = sol_EU.x(1);
    Beta_z = sol_EU.x(2);
    gamma_z = sol_EU.x(3);
    Beta_EU = sol_EU.x(4);
    gamma_EU = sol_EU.x(5);

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