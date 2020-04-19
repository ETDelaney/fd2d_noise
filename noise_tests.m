%%%
%You cannot generate band-limited "white" Gaussian noise. "White" noise means that the power spectral density is flat, which contradicts the notion of a passband.
%You can generate band-limited Gaussian noise. I'm not sure what you mean by "...signal ranging from 0 to 3 with a frequency of 0-6Hz", so I'll assume that you want a passband of 0 to 6 Hz. You did not tell us your sampling frequency, I'll assume 100 Hz.
%


Fs = 100;
d = fdesign.lowpass('Fp,Fst,Ap,Ast',6,10,0.5,40,Fs);
B = design(d);
% create white Gaussian noise the length of your signal
x = randn(1000,1);
% create the band-limited Gaussian noise
y = filter(B,x);
