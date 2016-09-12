
%==========================================================================
% compute sensitivity kernel for mu
%
% [K, u_adj_fft] = run3_adjoint( structure, noise_source, G_fft,
%       ref_station, rec, adjstf, wavefield_fwd, mode )
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% noise_source: contains spectrum and distribution of psd
% G_fft: Fourier transformed displacement Green function for ref. station
% ref_station: position of reference station
% rec: receiver positions, i.e. adjoint source locations
% adjstf: adjoint source time functions
% wavefield_fwd: forward wavefield
%
% mode: integer switch
%       == 0: do not Fourier transform adjoint wavefield
%       == 1: Fourier transform adjoint wavefield
%
% output:
%--------
% K: sensitivity kernels for source distribution and mu
% u_adj_fft: Fourier transformed adjoint wavefield - necessary for
%            structure kernel
%
%==========================================================================


function [K, u_adj_fft] = run3_adjoint(structure, noise_source, G_fft, ref_station, rec, adjstf, wavefield_fwd, mode)


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth, make_plots, plot_nth] = input_parameters();
    [X, Z, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);


    %- time axis ----------------------------------------------------------
    t = - (nt - 1) * dt:dt:(nt - 1) * dt;
    n_zero = nt;
    nt = length(t);


    %- prepare coefficients for Fourier transform -------------------------
    [~, n_sample, w_sample, dw, freq_samp] = input_interferometry();

    fft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    ifft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    for n = 1:nt
        for k = 1:n_sample
            fft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(- 1i * w_sample(k) * t(n)) * dt;
            ifft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(1i * w_sample(k) * t(n)) * dw;
        end
    end


    %- compute indices for receivers, i.e. the adjoint source locations ---
    n_receivers = size(rec, 1);
    rec_id = zeros(n_receivers, 2);

    for i = 1:n_receivers
        rec_id(i, 1) = find( min( abs(x - rec(i, 1)) ) == abs(x - rec(i, 1)), 1 );
        rec_id(i, 2) = find( min( abs(z - rec(i, 2)) ) == abs(z - rec(i, 2)), 1 );
    end


    %- initialise absorbing boundary taper a la Cerjan --------------------
    [absbound] = init_absbound();


    %- prepare figure for kernel building process -------------------------
    if( exist('OCTAVE_VERSION', 'builtin' ) == 0 )
        if (strcmp(make_plots, 'yes') && isempty(wavefield_fwd))
            
            fig = figure;
            set(fig, 'units', 'normalized', 'position', [0.13 0.11 0.60 0.75], 'Color', [1 1 1])
            cm = cbrewer('div', 'RdBu', 120);
            title_size = 14;
            font_size = 12;
            marker_size = 6;
            
            %- prepare forward wavefield plot -----------------------------
            ax1 = subplot(2, 2, 1);
            hold on
            xlabel(ax1, 'x [km]')
            ylabel(ax1, 'z [km]')
            colormap(ax1,cm)
            axis(ax1, 'image')
            box(ax1, 'on')
            set(ax1, 'LineWidth', 2, 'FontSize', font_size)
            
            
            %- prepare adjoint wavefield plot -----------------------------
            ax2 = subplot(2, 2, 2);
            hold on
            xlabel(ax2, 'x [km]')
            set(ax2, 'YTick', [])
            colormap(ax2,cm)
            axis(ax2,'image')
            box(ax2,'on')
            set(ax2, 'LineWidth', 2, 'FontSize', font_size)
            
            
            %- prepare kernel plot ----------------------------------------
            ax3 = subplot(2, 2, 3:4);
            hold on
            xlabel(ax3, 'x [km]')
            set(ax3, 'YTick', [])
            colormap(ax3,cm)
            axis(ax3,'image')
            box(ax3,'on')
            set(ax3, 'LineWidth', 2, 'FontSize', font_size)
            
            
            %- colorbar and title different for octave --------------------
            title(ax1, 'forward wavefield', 'FontSize', title_size)
            title(ax2, 'adjoint wavefield', 'FontSize', title_size)
            title(ax3, 'kernel build-up', 'FontSize', title_size)
            
            cb1 = colorbar('Position', [0.73, 0.11, 0.02, 0.34], ...
                'Ticks', [], 'AxisLocation', 'in', 'FontSize', 12);
            ylabel(cb1, 'fields and kernels are normalized', 'FontSize', font_size)
            colorbar('Position', [0.73, 0.11, 0.02, 0.34], ...
                'Ticks', [-1, 1], 'TickLabels', {'-', '+'}, 'AxisLocation', 'out', 'FontSize', font_size);
            
            
            %- get absorbing boundary width and initialize max-variables ------
            [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs();
            max_u = 0;
            max_M_tn = 0;
            max_interaction = 0;
            
            
            %- make movie -----------------------------------------------------
            % writerObj = VideoWriter('~/Desktop/test','MPEG-4');
            % writerObj.FrameRate = 6;
            % open(writerObj);
            
        end
        
    end


    %- allocate fields ----------------------------------------------------
    v = zeros(nx, nz);
    sxy = zeros(nx - 1, nz);
    szy = zeros(nx, nz - 1);
    u = zeros(nx, nz);

    K_mu = zeros(nx, nz);
    K_source = zeros(nx, nz);

    i_fw_in = 1;
    if (mode == 1)
        u_adj_fft = zeros(nx, nz, n_sample) + 1i * zeros(nx, nz, n_sample);
    else
        u_adj_fft = [];
    end


    %======================================================================
    % time loop
    %======================================================================

    %- only need first half of second adjoint run
    if (size(adjstf, 3) ~= 1)
        nt = n_zero;
    end


    for n = 1:nt


        %- compute divergence of current stress tensor --------------------
        DS = div_s(sxy, szy, dx, dz, nx, nz, order);


        %- add adjoint source time function -------------------------------
        if (size(adjstf, 3) == 1 && ~isempty(adjstf))

            for i = 1:n_receivers
                DS(rec_id(i, 1), rec_id(i, 2)) = DS(rec_id(i, 1), rec_id(i, 2)) + real(adjstf(i, n));
            end

        elseif (size(adjstf, 3) ~= 1 && ~isempty(adjstf))

            if (mod(n, freq_samp) == 0)
                T = zeros(nx, nz) + 1i * zeros(nx, nz);

                for k = 1:n_sample
                    T = T + noise_source.spectrum(k) * noise_source.distribution ...
                        .* conj(adjstf(:,:,k)) * ifft_coeff(n, k);
                end

                DS = DS + real(T);

            end

        end


        %- update velocity field ------------------------------------------
        v = v + dt * DS ./ structure.rho;


        %- apply absorbing boundary taper ---------------------------------
        v = v .* absbound;


        %- compute derivatives of current velocity and update stress tensor
        strain_dxv = dx_v(v, dx, nx, nz, order);
        strain_dzv = dz_v(v, dz, nx, nz, order);

        sxy = sxy + dt * structure.mu(1:nx-1, :) .* strain_dxv;
        szy = szy + dt * structure.mu(:, 1:nz-1) .* strain_dzv;


        %- calculate displacement and respective strain -------------------
        u = u + v * dt;
        strain_dxu = dx_v(u, dx, nx, nz, order);
        strain_dzu = dz_v(u, dz, nx, nz, order);


        %- build up source kernel -----------------------------------------
        if (~isempty(G_fft) && size(G_fft, 3) == n_sample && mod(n, freq_samp) == 0)

            M_tn = zeros(nx, nz) + 1i * zeros(nx, nz);

            for k = 1:n_sample
                M_tn = M_tn + noise_source.spectrum(k) .* G_fft(:,:,k) * ifft_coeff(n, k);
            end

            K_source = K_source + real(M_tn .* u);

        end


        %- build up structure kernel --------------------------------------
        if (~isempty(wavefield_fwd) && size(wavefield_fwd, 3) >= i_fw_in && mod(n, store_fwd_nth) == 0)

            K_mu(1:nx - 1,:) = K_mu(1:nx-1, :) - ...
                strain_dxu .* dx_v(wavefield_fwd(:,:, end - i_fw_in + 1), dx, nx, nz, order) * store_fwd_nth;

            K_mu(:, 1:nz - 1) = K_mu(:, 1:nz-1) - ...
                strain_dzu .* dz_v(wavefield_fwd(:,:, end - i_fw_in + 1), dz, nx, nz, order) * store_fwd_nth;

            i_fw_in = i_fw_in + 1;

        end


        %- save Fourier transformed adjoint state for second run ----------
        if (mode == 1)

            if (mod(n, freq_samp) == 0)
                for k = 1:n_sample
                    u_adj_fft(:,:, k) = u_adj_fft(:,:, k) + u * fft_coeff(n, k);
                end
            end

        end


        %- plot correlation wavefield -------------------------------------
        if( exist('OCTAVE_VERSION', 'builtin' ) == 0 )
            if (strcmp(make_plots, 'yes') && isempty(wavefield_fwd) ) % && n > 250)
                
                if (mod(n, plot_nth) == 0 || n == nt)
                    
                    %- kept plooting outside of function ----------------------
                    plot_kernel_build_up
                    drawnow
                    
                    %- make movie ---------------------------------------------
                    % M = getframe(fig);
                    % writeVideo(writerObj,M);
                     
                end
                
            end
        end


    end
    
    
    % if (strcmp(make_plots, 'yes'))
    %     close(writerObj);
    % end


    %- concatenate kernels ------------------------------------------------
    K = zeros(nx, nz, 2);
    K(:,:,1) = K_mu;
    K(:,:,2) = K_source;


end

