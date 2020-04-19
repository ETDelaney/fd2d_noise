
%==========================================================================
% compute Green function for reference station
% modified from original code: run1_forward_green so that can user can see
% the Green function plotted between two stations
%
% [ G_fft, G_out ] = plotgreen_run1_forward_green( structure, ref_station, mode )
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% ref_station: position of reference station - where we generate source for
% Green's function for first step in modelling correlation
% mode: integer switch
%       == 0: do not save wavefield
%       == 1: save wavefield
%
% output:
%--------
% G_fft: Fourier transformed displacement Green function for ref. station
% G_out: wavefield of displacement Green function
%
%==========================================================================


function [seismograms, G_fft, G_out] = plotgreen_run1_forward_green(structure, ref_station, rec, mode)

    %[G_fft, G_out] = run1_forward_green(structure, ref_station, mode)


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth, make_plots, plot_nth] = input_parameters();
    if strcmp(make_plots,'yes')
        [X, Z, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);
    else
        [~, ~, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);
    end
    
    % ###############3
    % delete this ------
    %- compute indices for receiver locations -----------------------------
    n_receivers = size(rec, 1);
    rec_id = zeros(n_receivers, 2);

    for i = 1:n_receivers
        rec_id(i, 1) = find( min( abs(x - rec(i, 1)) ) == abs(x - rec(i, 1)), 1 );
        rec_id(i, 2) = find( min( abs(z - rec(i, 2)) ) == abs(z - rec(i, 2)), 1 );
    end
    %##################3



    %- time axis ----------------------------------------------------------
    t = 0:dt:(nt - 1) * dt;


    %- prepare coefficients for Fourier transform -------------------------
    [~, n_sample, w_sample, ~, freq_samp] = input_interferometry();

    fft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    for n = 1:nt
        for k = 1:n_sample
            fft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(- 1i * w_sample(k) * t(n)) * dt;
        end
    end


    %- make source time function ------------------------------------------
   
    
    % test ricker
    dom_freq = 1;
    t0 = 1/(1.5*dom_freq);
    stdvar = 1/(2*pi*dom_freq);
    s=-(t-t0)/stdvar^2;
    s=s.*exp(-0.5*(t-t0).^2/stdvar^2);
    %s=-1500e9*ricker_wavelet(1/(2*pi*frequency),1/(1.5*frequency),t);
    stf = -1500*cumsum(s); % insert an integrated ricker
    
    % step function - to get G (since this code produces velocities)
    %stf = -1500*1.0e9 * ones(1, nt);

    
    
    % delta pulse - if we wanted dG/dt
    %stf = zeros(1, nt);
    %stf(1) = -1500*1.0e9;
    
    % dipole - if we wanted d2G/dt2
    %stf = zeros(1, nt);
    %stf(1) = 1.0e9;
    %stf(2) = -1*stf(1);
    
    % ricker - let us test this:
    %stf = zeros(1, nt);
    %delay = 2;
    %dom_freq = 1;
    %stf = (1-0.5*(2*pi*dom_freq)^2*(t-delay).^2);
    %stf = stf .* exp(-0.25*(2*pi*dom_freq)^2*(t-delay).^2);
    


    %- compute indices for source locations -------------------------------
    assert(size(ref_station, 1) == 1, ...
        'not possible to compute Green functions for more reference stations in one run of run1_forward_green.m')
    src_id = zeros(1, 2);
    src_id(1, 1) = find( min( abs(x - ref_station(1, 1)) ) == abs(x - ref_station(1, 1)), 1 );
    src_id(1, 2) = find( min( abs(z - ref_station(1, 2)) ) == abs(z - ref_station(1, 2)), 1 );


    %- initialise absorbing boundary taper a la Cerjan --------------------
    [absbound] = init_absbound();
    

    %- prepare figure for correlation wavefield ---------------------------
    if (strcmp(make_plots, 'yes'))

        fig = figure;
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5], 'Color', [1 1 1])
        title_size = 14;
        font_size = 12;
        marker_size = 6;
        
        ax1 = gca;
        hold(ax1, 'on')
        xlabel(ax1, 'x [km]')
        ylabel(ax1, 'z [km]')
        xlim(ax1, [0, Lx / 1000])
        ylim(ax1, [0, Lz / 1000])

        cm = cbrewer('div', 'RdBu', 120);
        colormap(cm)
        cb = colorbar;
        
        set(ax1, 'FontSize', font_size, 'position', [0.17, 0.204, 0.599, 0.624]);
        title(ax1, 'Green Function', 'FontSize', title_size)
        axis(ax1, 'image')
        box(ax1, 'on')
        set(ax1, 'LineWidth', 2)

        [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs();
        
        %- make movie -----------------------------------------------------
        % writerObj = VideoWriter('~/Desktop/test','MPEG-4');
        % writerObj.FrameRate = 6;
        % open(writerObj);

    end
    
    

    %- allocate fields ----------------------------------------------------
    v = zeros(nx, nz);
    sxy = zeros(nx - 1, nz);
    szy = zeros(nx, nz - 1);
    G_fft = zeros(nx, nz, n_sample) + 1i * zeros(nx, nz, n_sample);
    
    if strcmp(make_plots,'yes')
        u = zeros(nx, nz);
        seismograms = zeros(n_receivers, nt); % uncomment if you want to
        %see Green's function seismogram
    end

    i_fwd_out = 1;
    n_fwd = 0;
    for n = nt:(2 * nt - 1)
        if (mod(n, store_fwd_nth) == 0)
            n_fwd = n_fwd + 1;
        end
    end

    if (mode ~= 0)
        G_out = zeros(nx, nz, n_fwd, 'single');
    else
        G_out = single([]);
    end


    %======================================================================
    % time loop
    %======================================================================

    for n = 1:nt


        if (mode ~= 0 && mod(n + nt - 1, store_fwd_nth) == 0)
            G_out(:,:, i_fwd_out) = single(v);
            i_fwd_out = i_fwd_out + 1;
        end


        %- compute divergence of current stress tensor --------------------
        DS = div_s(sxy, szy, dx, dz, nx, nz, order);


        %- add point source -----------------------------------------------
        DS(src_id(1, 1), src_id(1, 2)) = DS(src_id(1, 1), src_id(1, 2)) + stf(n);


        %- update velocity field ------------------------------------------
        v = v + dt * DS ./ structure.rho;

        %- apply absorbing boundary taper ---------------------------------
        v = v .* absbound;


        %- compute derivatives of current velocity and update stress tensor
        strain_dxv = dx_v(v, dx, nx, nz, order);
        strain_dzv = dz_v(v, dz, nx, nz, order);

        sxy = sxy + dt * structure.mu(1:nx-1, :) .* strain_dxv;
        szy = szy + dt * structure.mu(:, 1:nz-1) .* strain_dzv;        
        

        %- accumulate Fourier transform of the displacement Greens function
        if (mod(n + nt - 1, freq_samp) == 0)

            for k = 1:n_sample
                G_fft(:,:,k) = G_fft(:,:,k) + v * fft_coeff(n, k);
            end

        end
        
        
        % ---- if we wish to plot the Green's function on-the-fly
        % ---- in time
        if (strcmp(make_plots, 'yes'))

        %- calculate displacement -----------------------------------------
        %    u = u + v * dt;
        % since we use a step function, you should not integrate the
        % velocity field, as you have techncially gotten your displacement
        % field ---> i.e., a delta pulse gives you dG/dt; where as a step
        % function gives you G when we use the first-order wave-propagation
        % setup. If you had a velocity recorder, it would record dG/dt for
        % a delta pulse.
        
        u = v; % we assume we are using a displacement recorder and that to
        % trick, we have given it a step function as the stf

        %- record seismograms ---------------------------------------------
        % uncomment if you want to
        % see Green's function seismogram
        for ir = 1:n_receivers
            seismograms(ir, n) = u(rec_id(ir, 1), rec_id(ir, 2));
        end


        %- plot correlation wavefield -------------------------------------
            if (n > 0.1 * nt && n < 0.9 * nt && mod(n, plot_nth) == 0)
            %if 1 == 1

                
                %- plot wavefield -----------------------------------------
                pcolor(ax1, X / 1000, Z / 1000, u')
                shading(ax1, 'interp')

                % we do not need to worry about source_type for Green's
                % Function calculation
                %if (strcmp(source_type, 'homogeneous'))
                %    m = 0.08;
                %elseif (strcmp(source_type, 'point'))
                %    m = 0.003;                    
                %elseif (strcmp(source_type, 'gaussian'))
                %    m = 0.15;
                %else
                    m = 0.4*max(max(abs(u)));
                %end
                caxis(ax1, [-m, m]);

                
                %- plot array ---------------------------------------------
                plot(ax1, ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ... 
                    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
                
                % we do not need this - we record Green's function for all
                % nodes: G(x, x_ref_station)
                plot(ax1,rec(:, 1) / 1000, rec(:, 2) / 1000, ....
                    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
                
                
                %- plot absorbing boundaries ------------------------------
                if(absorb_bottom); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [width, width] / 1000, 'k--'); end
                if(absorb_top); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [Lz - width, Lz - width] / 1000, 'k--'); end
                if(absorb_left); plot(ax1, [width, width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
                if(absorb_right); plot(ax1, [Lx - width, Lx - width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
                
                
                %- set labels ---------------------------------------------
                % if( exist('OCTAVE_VERSION', 'builtin' ) == 0 )
                %     set(cb, 'YTick', [-m m], 'TickLabels', {'-', '+'})
                % else
                    clabels = get(cb, 'YTick');
                    set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])
                % end
                
                xlabels = get(ax1, 'XTick');
                ylabels = get(ax1, 'YTick');
                set(ax1, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
                set(ax1, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
                xlim(ax1, [0, Lx / 1000])
                ylim(ax1, [0, Lz / 1000])
                
                
                % invoke plot ---------------------------------------------
                drawnow
                
                
                %- make movie ---------------------------------------------
                % M = getframe(fig);
                % writeVideo(writerObj,M);
                

            end

        end
        
        

    end
    test

end


