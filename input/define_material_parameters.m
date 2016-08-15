function [structure] = define_material_parameters( make_plots )

%==========================================================================
% generate material parameters mu [N/m^2] and rho [kg/m^3]
%==========================================================================
    

    % define make_plots if not specified
    if( nargin < 1 )
        make_plots = 'yes';
    end
    
    
    [Lx, Lz, nx, nz, ~, ~, ~, model_type] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    
    
    if( model_type == 1 )

        structure.rho = 3000.0 * ones(nx, nz);
        structure.mu = 4.8e10 * ones(nx, nz);        
    
    elseif( model_type == 2 )
        
        structure.rho = 3000.0 * ones(nx,nz);
        structure.mu = 4.8e10 * ones(nx,nz);
        
        x_sourcem = [2.0e5 2.3e5];
        z_sourcem = [2.0e5 2.0e5];
        x_width = [2.2e4 2.2e4];
        z_width = [2.2e4 5e4];
        strength = [-15.0e9, 8.0e9];
        
        for i = 1%:size( x_sourcem, 2 )
            structure.mu = structure.mu + strength(i) * exp( -( (X-x_sourcem(i)).^2 ) / x_width(i)^2 )' .* exp( -( (Z-z_sourcem(i)).^2 ) / z_width(i)^2 )' ;
        end
        
    end


    
    if( strcmp(make_plots, 'yes') )

        array = [];
        plot_models( sqrt(structure.mu./structure.rho), [], array, [0 0 0 0] );

    end

    
end


