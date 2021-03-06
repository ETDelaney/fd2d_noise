
function [Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters()

%==========================================================================
% set basic simulation parameters
%==========================================================================

% % tiny setup for gradient validation
Lx = 6.0e4;             % model extension in x-direction [m]
Lz = 6.0e4;             % model extension in y-direction [m]

nx = 50;                % grid points in x-direction
nz = 50;                % grid points in z-direction

dt = 0.09;              % time step [s]
nt = 150;                % number of iterations


% % small setup
% Lx = 4.0e5;         % model extension in x-direction [m]
% Lz = 4.0e5;         % model extension in y-direction [m]
% 
% nx = 300;           % grid points in x-direction
% nz = 300;           % grid points in z-direction
% 
% dt = 0.09;          % time step [s]
% nt = 900;           % number of iterations


% % setup Andreas
% Lx = 2.0e6;             % model extension in x-direction [m]
% % Lx = 3.0e6;             % model extension in x-direction [m]
% Lz = 2.0e6;             % model extension in y-direction [m]
% 
% nx = 600;               % grid points in x-direction
% % nx = 900;               % grid points in x-direction
% nz = 600;               % grid points in z-direction
% 
% dt = 0.23;              % time step [s]
% nt = 1300;              % number of iterations
% % nt = 1600;              % number of iterations
% % nt = 3000;            % number of iterations


order = 4;              % finite-difference order (2 or 4)


%==========================================================================
% set basic inversion parameters
%==========================================================================

% store every nth time step of forward field
fw_nth = 1;


%==========================================================================
% model type
%==========================================================================

model_type = 1;
% model_type = 200;
% model_type = 666;
% model_type = 888;
% model_type = 999;


% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=different layered medium
% 6=vertical gradient medium in mu
% 7=another layered medium
% 666 = put source picture in ../models-folder
% 999 = IUGG structure, two pills


%==========================================================================
% source type
%==========================================================================

source_type = 'homogeneous';
% source_type = 'gaussian';
% source_type = 'point';

% number of frequency bands
n_basis_fct = 0;  

