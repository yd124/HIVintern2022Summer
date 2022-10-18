data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);

tic
out = simulate(1000,dpi,y,0.1*num);



tb = array2table(out,...
            'VariableNames', ...
            {'J','b','dlt','p','d'});
writetable(tb, '2sim100v1.csv');

toc




%% find the rough minimum error 
function out = simulate(max_N,dpi,y,erbd)
    
    h = 0.01;
    ti = 0:h:dpi(end);

    out = zeros(max_N,5);
    N = 0;
    n = 0;

    while N < max_N
        n = n + 1;
        fprintf('n %d\n',n);

        %% Unif(a,b) -> a+(b-a)*rand
        %% N(mu,sigma) -> mu+sigma*randn
        init = [10^4 0 10^-3];
        
        b = 10^-6*rand;     

        p = 10^5*rand;      
        
%         init = [10^4 0 10^-3];
        
        dlt = 1*rand;
        d = 1*rand; 

%         b = 10^-6 + (10^-5-10^-6)*rand;
%         p = 10^2.6 + (10)

        if b <= 0 || dlt <= 0 || p <= 0 || d <= 0
            continue
        end 
        

        try
            xa = pred(ti,init,b,dlt,p,d);
        catch
            fprintf('signluar points');
            continue
        end
         
        V = xa(:,3);
        if min(V) <= 0
            continue
        end 
      
        y_hat = log10(V(dpi/h+1));
        error = MSE(y,y_hat);
        if error < erbd
            N = N + 1;
            fprintf('N %d\n',N);
            out(N,1) = error;
            out(N,2:end) = log10([b dlt p d]);
        end 
    end 
end 






function out = pred(ti,init,b,dlt,p,d)

    f = @(t,x) [d*(init(1)-x(1))-b*x(1)*x(3);...
            b*x(1)*x(3)-dlt*x(2);...
            p*x(2)-23*x(3) ];    
    [t,xa] = ode45(f,ti,init);
    out = xa;

end




%% Mean-squared error
function out = MSE(y,y_hat)
    out = 1/size(y,1) *sum((y-y_hat).^2);
end
