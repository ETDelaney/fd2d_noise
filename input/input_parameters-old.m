
%==========================================================================
% set simulation parameters, model and source type, and plotting behaviour
%
% (see function body for more explanation)
%
%==========================================================================


function [Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, store_fwd_nth, make_plots, plot_nth] = input_parameters()


    %======================================================================
    % basic simulation parameters
    %======================================================================

%     Lx = 4.0e5;             % model extension in x-direction [m]
%     Lz = 4.0e5;             % model extension in y-direction [m]
%     nx = 300;               % grid points in x-direction
%     nz = 300;               % grid points in z-direction
% 
%     dt = 0.09;              % time step [s]
%     nt = 800;               % number of time steps
% 
% 
%     order = 4;              % finite-difference order (2 or 4)

    Lx = 4e3;             % model extension in x-direction [m]
    Lz = 4.5e3;             % model extension in y-direction [m]
    nx = 125;               % grid points in x-direction
    nz = 141;               % grid points in z-direction

    dt = 0.025;              % time step [s]
    nt = 600;  %500             % number of time steps

    order = 4;              % finite-difference order (2 or 4)



    %======================================================================
    % model type
    %======================================================================

    model_type = 100;
    % 1 = homogeneous 
    % 2 = slow anomaly in the center of the domain
    % 100 = starting model for SWIM array for freq. 0.55-0.95 Hz
    % 101 = test perturbation


    %======================================================================
    % source type
    %======================================================================

    source_type = 'homogeneous';
    %source_type = 'point';
    %source_type = 'gaussian';


    %======================================================================
    % plotting parameters
    %======================================================================

    make_plots = 'no';      % 'yes' or 'no'
    plot_nth = 10;          % plot every nth time step
    
    
    %======================================================================
    % parameter for structure kernel calculation
    %======================================================================

    store_fwd_nth = 4;      % store forward wavefield every nth time step


end

