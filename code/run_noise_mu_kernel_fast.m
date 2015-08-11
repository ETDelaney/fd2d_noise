function [X,Z,K_mu_final] = run_noise_mu_kernel_fast(C_2_dxv, C_2_dzv, mu, stf, adsrc)

%==========================================================================
% run simulation to compute sensitivity kernel for rho and mu
% (only one-sided)
% fast means ready for conversion to mex-files
%
% input:
%--------
% C_2_dxv & C_2_dzv: strain of correlation wavefield
% mu [N/m^2]
% stf: adjoint source time function
% adsrc: adjoint source positions
%
% output:
%--------
% X, Z: coordinate axes
% K_mu_final: sensitivity kernel for mu
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
nt = length(t);


%- reverse adjoint source time function -----------------------------------
stf = fliplr(stf);


%- compute indices for adjoint source locations ---------------------------    
ns = size(adsrc,1);
adsrc_id = zeros(ns,2);
for i=1:ns
    adsrc_id(i,1) = min( find( min(abs(x-adsrc(i,1))) == abs(x-adsrc(i,1))) );
    adsrc_id(i,2) = min( find( min(abs(z-adsrc(i,2))) == abs(z-adsrc(i,2))) );
end


%- initialise interferometry ----------------------------------------------        
[~,n_sample,w_sample,dw] = input_interferometry();

G_1_dxv = zeros(nx-1,nz,n_sample) + 1i*zeros(nx-1,nz,n_sample);
G_1_dzv = zeros(nx,nz-1,n_sample) + 1i*zeros(nx,nz-1,n_sample);
            

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
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    for i=1:ns
        DS( adsrc_id(i,1),adsrc_id(i,2) ) = DS( adsrc_id(i,1),adsrc_id(i,2) ) + stf(i,n);
    end
    
    
    %- update velocity field ----------------------------------------------    
    v = v + dt*DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---
    strain_dxv = dx_v(v,dx,dz,nx,nz,order);
    strain_dzv = dz_v(v,dx,dz,nx,nz,order);
    
    sxy = sxy + dt*mu(1:nx-1,:) .* strain_dxv;
    szy = szy + dt*mu(:,1:nz-1) .* strain_dzv;
     
    
    %- accumulate Fourier transform of the velocity field -----------------
    if( mod(n,5) == 1 )
        for k=1:n_sample           
            G_1_dxv(:,:,k) = G_1_dxv(:,:,k) + strain_dxv(:,:) * exp(-1i*w_sample(k)*t(n)) * dt;
            G_1_dzv(:,:,k) = G_1_dzv(:,:,k) + strain_dzv(:,:) * exp(-1i*w_sample(k)*t(n)) * dt;
        end
    end

end


%==========================================================================
% compute structure kernels 
%==========================================================================

%- accumulate kernel by looping over frequency
K_mu = zeros(nx,nz) + 1i*zeros(nx,nz);
for k=1:n_sample   
    K_mu(1:nx-1,:) = K_mu(1:nx-1,:) + G_1_dxv(:,:,k) .* C_2_dxv(:,:,k) / w_sample(k)^2 * dw;
    K_mu(:,1:nz-1) = K_mu(:,1:nz-1) + G_1_dzv(:,:,k) .* C_2_dzv(:,:,k) / w_sample(k)^2 * dw;
end
K_mu_final = real(K_mu);

