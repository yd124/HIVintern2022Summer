
cd("/Users/yuqingdai/Documents/GitHub/HIVintern2022Summer/Old Paper/MATLAB/profile_likelihood")

LSEs = table2array(readtable('LSE1.csv'));
sampled_params = table2array(readtable('sampled_params.csv'));
% original_params = table2array(readtable('original_params.csv'));
results = table2array(readtable('../new_work/results.csv'));
original_params = table2array(readtable('original_params.csv'));


xx_error = original_params(1);
xx = original_params(2:end);

num_params = size(LSEs,1);

names = ["b0" "bi" "k" "dlt" "p" "d" "tau"];

tiledlayout(4,2);
for ip = 1:num_params
    nexttile
    p1 = sampled_params(ip,:);
    a1 = LSEs(ip,:);
%     scatter(p1(a1>0),a1(a1>0),'filled');
    plot(p1(a1>0),a1(a1>0));
    xlabel(names(ip));
    ylabel("LSE");
    hold on
    for i = 1:size(results,1)
        scatter(results(i,1+ip),results(i,1),'filled');
        hold on
    end
    legend('','x1','x2','x3', 'x4','x5','x6','Location','northwest')

    hold off
end 