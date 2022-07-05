data = readtable('sim100v5.csv');
y = data.J;

c = ismember(y,mink(y,10));
idx = find(c);
disp(data(idx,:))


writetable(data(idx,:), 'first10.csv');