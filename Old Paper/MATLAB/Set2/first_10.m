data = readtable('sim100v2.csv');
y = data.J;

c = ismember(y,mink(y,10));
idx = find(c);
disp(data(idx,:))


writetable(data(idx,:), 'first10v2.csv');