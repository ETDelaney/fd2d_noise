
addpath(genpath('../../'))

% G_2, C_2, C_2_dxv, C_2_dzv
i1 = coder.typeof(complex(zeros(2,2,2)),[inf,inf,inf],1);

% mu, adstf
i2 = coder.typeof(zeros(2,2),[inf,inf],1);

% src, rec
i3 = coder.typeof(zeros(2,2),[inf,2],1);

% integer switch
i4 = coder.typeof(1,1,0);

% source dist
i5 = coder.typeof(zeros(2,2,2),[inf,inf,inf],1);

% run_forward_green_fast(mu,src)
codegen run_forward_green_fast.m -args {i2,i3}

% run_forward_correlation_fast(G_2, source_dist, mu, rec)
codegen run_forward_correlation_fast.m -args {i1,i5,i2,i3,i4}


rmpath(genpath('../../'))
clear i1
clear i2
clear i3
clear i4
clear i5