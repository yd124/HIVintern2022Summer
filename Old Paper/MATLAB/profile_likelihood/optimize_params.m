cd('~/Documents/GitHub/HIVintern2022Summer/Old Paper/MATLAB/profile_likelihood/')
% @pred.m; 


data = readtable('../../data/data1/logvRNA.csv');
patients = unique(data.patient,'stable');

x1 = [0.409E-6, 0.233E-6, 0.249, 0.775, 14.5E3, 0.03, 7];
x2 = [0.431E-6, 0.140E-6, 0.077, 0.420, 10E3, 0.021, 24];
x3 = [0.201E-6, 0.001E-6, 0.013, 1.048, 30.172E3, 0.036, 10];
x4 = [9.203E-6, 0.011E-6, 0.013, 0.851, 0.548E3, 0.055, 12];
x5 = [0.485E-6, 0.291E-6, 0.096, 0.803, 11.425E3, 0.033, 5];
x6 = [0.057E-6, 0.004E-6, 0.021, 0.821, 89.892E3, 0.003, 22];

x = [x1;x2;x3;x4;x5;x6];

% parfor loop can use cpu and speed it up
for i = 1:length(patients)
% parfor i = 1:length(patients)
    xx = x(i,:);
    sub_data = data(data.patient == string(patients(i)),:);
    dpi = sub_data.dpi;
    
    y = sub_data.log_vRNA;
    y_hat = pred(xx,dpi);
    disp(y_hat)
%     [opt_params,fval] = best_param(xx,dpi,y);
%     disp(fval);

end
