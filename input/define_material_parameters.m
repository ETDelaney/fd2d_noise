
%==========================================================================
% generate material parameters mu [N/m^2] and rho [kg/m^3]
%
% [structure] = define_material_parameters( make_plots )
%
% input:
%--------
% make_plots: 'yes' or 'no' (optional)
%
% output:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
%
%==========================================================================


function [structure] = define_material_parameters(make_plots)


    %- define make_plots if not specified ---------------------------------
    if (nargin < 1)
        make_plots = 'yes';
    end
    
    
    %- check path ---------------------------------------------------------
    fd2d_path();


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, ~, ~, ~, model_type] = input_parameters();
    [X, Z] = define_computational_domain(Lx, Lz, nx, nz);


    %- set up different models --------------------------------------------
    if (model_type == 1)

        structure.rho = 3000.0 * ones(nx, nz);
        structure.mu = 4.8e10 * ones(nx, nz);

        
    elseif (model_type == 2)

        structure.rho = 3000.0 * ones(nx, nz);
        structure.mu = 4.8e10 * ones(nx, nz);

        x_sourcem = [2.0e5, 2.3e5];
        z_sourcem = [2.0e5, 2.0e5];
        x_width = [2.2e4, 2.2e4];
        z_width = [2.2e4, 5e4];
        strength = [- 15.0e9, 8.0e9];

        for i = 1 %:size( x_sourcem, 2 )
            structure.mu = structure.mu + strength(i) * exp(- ((X - x_sourcem(i)) .^ 2) / x_width(i) ^ 2)' .* exp(- ((Z - z_sourcem(i)) .^ 2) / z_width(i) ^ 2)';
        end


        %==================================================================
        % define your own velocity model, e.g.
        %
        % elseif( model_type == 3 )
        %
        %     structure.rho = 3000.0 * ones(nx, nz) + ...;
        %     structure.mu = 4.8e10 * ones(nx, nz) + ...;
        %
        %==================================================================
        
        
    elseif (model_type == 100)

        % initial guess for Scholte Waves in frequency range 0.55-0.95 Hz
        structure.rho = 1800.0 * ones(nx, nz); % based on shallow analysis at OL and Statfjord; should check Oseberg
        structure.mu = 2.33e8 * ones(nx, nz); % adjust to make velocity approximate 360 m/s 
    
    elseif (model_type == 101)

        % initial guess for Scholte Waves in frequency range 0.55-0.95 Hz
        structure.rho = 1800.0 * ones(nx, nz); % based on shallow analysis at OL and Statfjord; should check Oseberg
        structure.mu = 2.33e8 * ones(nx, nz); % adjust to make velocity approximate 360 m/s 
        
        % now perturb from start model
        scale_factor = 15000000;
        load('K_0.55-0.95_100-0.mat','gradient_masked')
        %structure.rho = structure.rho + scale_factor*gradient_masked(:,:,1);
        
% original tests using synthetics        
%     elseif (model_type == 101)
% 
%         % initial guess for Scholte Waves in frequency range 0.55-0.95 Hz
%         structure.rho = 1800.0 * ones(nx, nz); % based on shallow analysis at OL and Statfjord; should check Oseberg
%         structure.mu = 2.33e8 * ones(nx, nz); % adjust to make velocity approximate 360 m/s 
%         
%         x_sourcem = [2.0e3, 2.3e3];
%         z_sourcem = [2.0e3, 2.0e3];
%         x_width = [2.2e2, 2.2e2];
%         z_width = [2.2e2, 5e2];
%         strength = [1500, 8.0e7];
% 
%         for i = 1 %:size( x_sourcem, 2 )
%             %structure.mu = structure.mu + strength(i) * exp(- ((X - x_sourcem(i)) .^ 2) / x_width(i) ^ 2)' .* exp(- ((Z - z_sourcem(i)) .^ 2) / z_width(i) ^ 2)';
%             structure.rho = structure.rho + strength(i) * exp(- ((X - x_sourcem(i)) .^ 2) / x_width(i) ^ 2)' .* exp(- ((Z - z_sourcem(i)) .^ 2) / z_width(i) ^ 2)';
%         end
        
        
    end


    %- plot model ---------------------------------------------------------
    if (strcmp(make_plots, 'yes'))

        plot_models(sqrt(structure.mu ./ structure.rho), [], [], [0, 0, 0, 0]);

    end


end


