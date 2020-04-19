% taking the raw data, shape into the process that allows us to make it
% work with the f2d2 code - we have to make redundant pairs (already taken
% care of in the python code) - but then reshaped from 29412 to 172,172...
% with the samples according to the length of the time vector of the
% modelling code

load('C:\Users\Tiberius\Dropbox\PhD\fd2d_noise_05jan2020\output\2014.05.month_0.55.to.0.95_3_redundant.mat')
figure
plot(cc_matrix(1,:))
cc_data = flip(cc_matrix);
figure
plot(cc_data(1,:))

for i = 1:29412
cc_data(i,:) = flip(cc_matrix(i,:));
end
figure
plot(cc_data(1,:))
cc_data_resample2x = resample(cc_data',2,1)';
figure
plot(cc_data_resample2x(1,:))

cc_data_final = cc_data_resample2x(:,[602:1800]);
cc_fin_final = zeros(172,171,1199);
for i=1:172
cc_fin_final(i,:,:) = cc_data_final(1+171*(i-1):171*i,:);
end
correlations = cc_fin_final;