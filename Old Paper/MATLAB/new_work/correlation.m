k_tb = table2array(readtable('k.csv'));

% ab_tb = readtable('../../data/data1/antibodies_rate.csv');
% writetable(ab_tb,'ab.csv');

ab_tb = readtable('ab.csv');

subplot(1,3,1)
scatter(ab_tb.IgM,k_tb(:,1),'filled')
xlabel('Slope of IgM increase');
ylabel('Slope of infectivity decay');


subplot(1,3,2)
scatter(ab_tb.IgG,k_tb(:,1),'filled')
xlabel('Slope of IgG increase');
ylabel('Slope of infectivity decay');


subplot(1,3,3)
scatter(ab_tb.IgM_IgG,k_tb(:,1),'filled')
xlabel('Slope of IgM+IgG increase');
ylabel('Slope of infectivity decay');
