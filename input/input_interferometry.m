
%==========================================================================
% input parameters for Fourier transform
%
% [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()
%
% output:
%--------
% f_sample: frequency vector [Hz]
% n_sample: number of frequency samples
% w_sample: frequency vector converted to angular frequency
% dw: step in angular frequency
% freq_samp: calculate Fourier transform every freq_samp time step
%
%==========================================================================


function [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()


    %f_sample = 0.02:0.004:0.2;      % in Hz
    %f_sample = 0.3:0.01:1.2; % this is a bit broad, no?
    %f_sample = 0.35:0.01:1.15;
    f_sample = 0:0.0166527:1.665; % experiment with a forced spectrum
    %f_sample = 0:0.014275517487508882:1.67023555;
    
    freq_samp = 1; %5;


    %- derived quantities -------------------------------------------------
    n_sample = length(f_sample);
    w_sample = 2 * pi * f_sample;
    dw = w_sample(2) - w_sample(1);


end
