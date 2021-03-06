data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);

tic
% out = simulate(10,num,dpi,y,Inf);
% erbd = min(out(:,1));


out = simulate(1000,dpi,y,0.1*num);

% y = out(:,1);
% 
% c = ismember(y,mink(y,1));
% idx = find(c);
% disp(out(idx,:))



tb = array2table(out,...
            'VariableNames', ...
            {'J','b0','bi','k','dlt','p','d','tau'});
writetable(tb, 'sim100v1.csv');

toc




%% find the rough minimum error 
function out = simulate(max_N,dpi,y,erbd)
    
    h = 0.01;
    ti = 0:h:dpi(end);

    out = zeros(max_N,8);
    N = 0;
    n = 0;

    while N < max_N
        n = n + 1;
%         fprintf('n %d\n',n);

        %% Unif(a,b) -> a+(b-a)*rand
        %% N(mu,sigma) -> mu+sigma*randn
    
        
        bi = 10^-5*rand;     
        b0 = bi + 10^-5*rand;
        
%         b0 = 10^-4*rand;
        if b0 <= bi
            continue 
        end

        p = 10^4*rand;      
        
        init = [10^4 0 10^-3];
        
        dlt = min(5,init(1)*p*bi/23)*rand;
        k = 30*rand;   
        d = 1*rand; 

        if bi <= 0 || k <= 0 || dlt <= 0 || p <= 0 || d <= 0
            continue
        end 
        
        tau = 30*rand;
        try
            xa = pred(ti,init,b0,bi,k,dlt,p,d,tau);
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
            out(N,2:end) = log10([b0 bi k dlt p d tau]);
        end 
    end 
end 




% if N > 0
%     plot(ti, best_logV);
%     hold on
%     scatter(dpi,y);
%     hold off
% end 


% params = [0.409*10^-6, 0.233*10^-6, 0.249, 0.775, 14.5*10^3, 0.03, 7];
% logV = pred_logV(params, ts, init);



%% decay function 
function out=b(t,b0,bi,k,tau)
    if t <= tau
        out = b0;
    else
        out = bi+(b0-bi)*exp(-k*(t-tau));
    end
end

function out = pred(ti,init,b0,bi,k,dlt,p,d,tau)

%     b0 = params(1);
%     bi = params(2);
%     k = params(3);
%     dlt = params(4);
%     p = params(5);
%     d = params(6);
%     tau = params(7);

    f = @(t,x) [d*(init(1)-x(1))-b(t,b0,bi,k,tau)*x(1)*x(3);...
            b(t,b0,bi,k,tau)*x(1)*x(3)-dlt*x(2);...
            p*x(2)-23*x(3) ];    
%     options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [t,xa] = ode45(f,ti,init);
    out = xa;

end




%% Mean-squared error
function out = MSE(y,y_hat)
    out = 1/size(y,1) *sum((y-y_hat).^2);
end
