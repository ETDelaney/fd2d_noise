function [X,Z,K_s] = run_noise_source_kernel_fast(G_2,mu,stf,adsrc)

%==========================================================================
% run simulation to compute sensitivity kernel for noise power-spectral
% density distribution
% fast means ready for conversion to mex-files
%
% input:
%--------
% G_2: Green function of reference station
% mu [N/m^2]
% stf: adjoint source time function
% adsrc: adjoint source positions
%
% output:
%--------
% X, Z: coordinate axes
% K_s: sensitivity kernel
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[~,rho] = define_material_parameters(nx,nz,model_type); 
mu = reshape(mu,nx,nz);

%- time axis --------------------------------------------------------------    
t = -(nt-1)*dt:dt:(nt-1)*dt;
    

%- compute indices for adjoint source locations ---------------------------
ns = size(adsrc,1);
adsrc_id = zeros(ns,2);

for i=1:ns
    adsrc_id(i,1) = min(find(min(abs(x-adsrc(i,1)))==abs(x-adsrc(i,1))));
    adsrc_id(i,2) = min(find(min(abs(z-adsrc(i,2)))==abs(z-adsrc(i,2))));
end


%- initialise interferometry ----------------------------------------------       
[~,n_sample,w_sample] = input_interferometry();

G_1 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
K_s = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
           

%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

for n=1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    for i=1:ns
        DS(adsrc_id(i,1),adsrc_id(i,2)) = DS(adsrc_id(i,1),adsrc_id(i,2)) + stf(i,n);
    end
    
    
    %- update velocity field ----------------------------------------------    
    v = v + dt*DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---    
    sxy = sxy + dt*mu(1:nx-1,:) .* dx_v(v,dx,dz,nx,nz,order);
    szy = szy + dt*mu(:,1:nz-1) .* dz_v(v,dx,dz,nx,nz,order);
     
    
    %- accumulate Fourier transform of the velocity field -----------------
    if( mod(n,5) == 1 )
        for k=1:n_sample
            G_1(:,:,k) = G_1(:,:,k) + v(:,:) * exp(-1i*w_sample(k)*t(n))*dt;
        end
    end    

end


%==========================================================================
% compute noise source kernels as function of frequency 
%==========================================================================

%- multiply Fourier transformed Greens functions
for k=1:n_sample
    K_s(:,:,k) = G_1(:,:,k) .* conj(G_2(:,:,k)) / (1i*w_sample(k)+eps);
end

K_s = real(K_s);

