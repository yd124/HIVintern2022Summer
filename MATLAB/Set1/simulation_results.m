data = readtable('sim100v6.csv');
y = data.J;


idx = find(y==min(y));
disp(data(idx,:))

% tb = array2table(data(idx,:),...
%             'VariableNames', ...
% %             {'J','b0','bi','k','dlt','p','d','tau'});
% writetable(tb, 'params_sim100v3.csv');

tiledlayout(3,3);

nexttile
scatter(data.b0,y,'filled')
xlabel('b0');
ylabel('J');

nexttile
scatter(data.bi,y,'filled')
xlabel('bi');
ylabel('J');

nexttile
scatter(data.k,y,'filled')
xlabel('k');
ylabel('J');

nexttile
scatter(data.dlt,y,'filled')
xlabel('dlt');
ylabel('J');


nexttile
scatter(data.p,y,'filled')
xlabel('p');

nexttile
scatter(data.d,y,'filled')
xlabel('d');
ylabel('J');


nexttile
scatter(data.tau,y,'filled')
xlabel('tau');
ylabel('J');

