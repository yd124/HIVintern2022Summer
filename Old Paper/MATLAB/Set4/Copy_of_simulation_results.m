data = readtable('sim1000v2.csv');
y = data.J;


idx = find(y==min(y));
disp(data(idx,:))

tiledlayout(2,2);

nexttile
scatter(data.b,y,'filled')
xlabel('b');
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


