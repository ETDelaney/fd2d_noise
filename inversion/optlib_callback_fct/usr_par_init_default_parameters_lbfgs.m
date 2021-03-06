
function [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par)
    
    [~,~,nx,nz,~,~,~,~,~,n_basis_fct,~] = input_parameters();
    if( isfield( usr_par, 'config') )
        if( ~isfield( usr_par.config, 'nx') )
            usr_par.config.nx = nx;
        end
        if( ~isfield( usr_par.config, 'nz') )
            usr_par.config.nz = nz;
        end
        if( ~isfield( usr_par.config, 'n_basis_fct') )
            usr_par.config.n_basis_fct = n_basis_fct;
        end
    else
        usr_par.config.nx = nx;
        usr_par.config.nz = nz;
        usr_par.config.n_basis_fct = n_basis_fct;
    end
    

    if( ~isfield( usr_par, 'cluster') )
        usr_par.cluster = 'local';
    end

    
    if( ~isfield( usr_par, 'type') )
        usr_par.type = 'joint';
    end
    
    
    if( ~isfield( usr_par, 'use_mex') )
        usr_par.use_mex = 'no';
    end
    
    
    if( ~isfield( usr_par, 'verbose') )
        usr_par.verbose = 'no';
    end
    
    if( strcmp( usr_par.use_mex, 'no' ) )
        ! cp ../code/run_forward1_green.m ../code/mex_functions/run_forward1_green_mex.m
        ! cp ../code/run_forward2_correlation.m ../code/mex_functions/run_forward2_correlation_mex.m
        ! cp ../code/run_noise_adjoint.m ../code/mex_functions/run_noise_adjoint_mex.m
    end
    
    
    if( isfield( usr_par, 'initial') )
        if( ~isfield( usr_par.initial, 'ref_source') )
            usr_par.initial.ref_source = 0;
        end
        if( ~isfield( usr_par.initial, 'ref_structure') )
            usr_par.initial.ref_structure = 1;
        end
        if( ~isfield( usr_par.initial, 'mu_0') )
            usr_par.initial.mu_0 = 4.8e10;
        end
    else
        usr_par.initial.ref_source = 0;
        usr_par.initial.ref_structure = 1;
        usr_par.initial.mu_0 = 4.8e10;
    end

    
    
    if( isfield( usr_par, 'measurement') )
        if( ~isfield( usr_par.measurement, 'source') )
            usr_par.measurement.source = 'waveform_difference';
        end
        if( ~isfield( usr_par.measurement, 'structure') )
            usr_par.measurement.structure = 'waveform_difference';
        end
    else
        usr_par.measurement.source = 'waveform_difference';
        usr_par.measurement.structure = 'waveform_difference';
    end

    
    if( isfield(usr_par,'kernel') )
        % if( isfield( usr_par.kernel, 'imfilter') )
        %     if( ~isfield( usr_par.kernel.imfilter, 'source' ) )
        %         usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
        %     end
        %     if( ~isfield( usr_par.kernel.imfilter, 'structure' ) )
        %         usr_par.kernel.imfilter.structure = fspecial('gaussian',[1 1], 1);
        %     end            
        % else 
        %     usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
        %     usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;
        % end
        
        
        if( isfield( usr_par.kernel, 'sigma') )
            if( ~isfield( usr_par.kernel.sigma, 'source' ) )
                usr_par.kernel.sigma.source = [1 1];
            end
            if( ~isfield( usr_par.kernel.sigma, 'structure' ) )
                usr_par.kernel.sigma.structure = [1 1];
            end
        else
            usr_par.kernel.sigma.source = [1 1];
            usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;
        end
        
        
        if( ~isfield( usr_par.kernel, 'weighting') )
            usr_par.kernel.weighting = 0.5;
        end
    else
        % usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
        % usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;
        
        usr_par.kernel.sigma.source = [1 1];
        usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;
        
        usr_par.kernel.weighting = 0.5;
    end
    
    
    if( isfield( usr_par, 'ring' ) )
        if( ~isfield( usr_par.ring, 'switch' ) )
            usr_par.ring.switch = 'no';
        end
    else
        usr_par.ring.switch = 'no';
    end
    
    
    if( ~isfield( usr_par, 'network') )
        error('\nload array!\n')
    end
    
    
    if( ~isfield( usr_par, 'data') )
        error('\nload data!\n')
    end
    
    
    if( ~isfield( usr_par, 'veldis') )
        usr_par.veldis = 'dis';
    end

    
    if( isfield(usr_par,'filter') )
        if( ~isfield( usr_par.filter, 'apply_filter') )
            usr_par.filter.apply_filter = 'no';
        end
    else
        usr_par.filter.apply_filter = 'no';
    end


    if( isfield(usr_par,'regularization') )
        if( ~isfield( usr_par.regularization, 'alpha') )
            usr_par.regularization.alpha = 0.0;
        end
        
        if( ~isfield( usr_par.regularization, 'beta') )
            usr_par.regularization.beta = 0.0;
        end
        
        if( ~isfield( usr_par.regularization, 'weighting') )
            usr_par.regularization.weighting = 1.0;
        end
    else
        usr_par.regularization.alpha = 0.0;
        usr_par.regularization.beta = 0.0;
        usr_par.regularization.weighting = 1.0;
    end


end
