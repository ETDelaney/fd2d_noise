
addpath(genpath('../../'))


% G_2, adstf
i1 = coder.typeof(complex(zeros(2,2,2)), [inf,inf,inf], 1);

% mu, rho, spectrum, dmu
i2 = coder.typeof(zeros(2,2), [inf,inf], 1);

% src, rec
i3 = coder.typeof(zeros(2,2), [inf,2], 1);

% integer switch
i4 = coder.typeof(1, 1, 0);

% source_dist
i5 = coder.typeof(zeros(2,2,2), [inf,inf,inf], 1);

% u_fwd, G_in, C_in, u_in
i6 = coder.typeof((zeros(2,2,2,'single')), [inf,inf,inf], 1);



% run_forward1_green( mu, rho, src, rec, mode, dmu, G_in )
codegen run_forward1_green.m -args {i2, i2, i3, i3, i4, i2, i6}


% run_forward2_correlation( mu, rho, G_fft, spectrum, source_dist, rec, mode, dmu, C_in )
codegen run_forward2_correlation.m -args {i2, i2, i1, i2, i5, i3, i4, i2, i6}


% run_noise_adjoint( mu, rho, u_fwd, adstf, adsrc, rec, spectrum, source_dist, G_fft, mode, dmu, u_in )
codegen run_noise_adjoint.m -args {i2, i2, i6, i1, i3, i3, i2, i5, i1, i4, i2, i6}



rmpath(genpath('../../'))