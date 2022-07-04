data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);

tic
% out = simulate(10,num,dpi,y,Inf);
% erbd = min(out(:,1));


out = simulate(100,num,dpi,y,0.1);

tb = array2table(out,...
            'VariableNames', ...
            {'J','b0','bi','k','dlt','p','d','tau'});
writetable(tb, 'sim100v2.csv');

toc




%% find the rough minimum error 
function out = simulate(max_N,num,dpi,y,erbd)
    
    h = 0.01;
    out = zeros(max_N,8);
    N = 0;

    n = 0;
    while N < max_N
        n = n + 1;
        fprintf('n %d\n',n);
        
        

        %% Unif(a,b) -> a+(b-a)*rand
        bi = 10^-6*rand;
        b0 = bi + 10^-6*rand;

        
        p = 10^2 + (10^5-10^2)*rand; 
      
        init = [10^4 0 10^-3];

        dlt = init(1)*p*bi/23*rand;

        k = 20*rand;   
        d = 2*rand; 
        
        if bi == 0 || k == 0 || dlt == 0 || p == 0 || d == 0
            continue
        end 
        

        tau = unidrnd(30);
%         params = [b0 bi k dlt p d tau];
        
    
        start = 0;
        V = zeros(dpi(end)/h+1,1);
     
        for i = 1:num
    
            sti = start:h:dpi(i);
            xa = pred(sti,init,b0,bi,k,dlt,p,d,tau);
    
            Vi = xa(end,3);
    
    
            if (Vi > 0) && (abs(log10(Vi)-y(i))< 1)
    
                try
                    idx = round(sti/h+1);
                    V(idx) = xa(:,3);
                catch
                    fprintf('signluar points');
                    break;
                end 
    
                start = dpi(i);
                init = xa(end,:);
    
            else 
                break
            end
    
        end 
        
        if i == num 
            y_hat = log10(V(dpi/h+1));
            error = MSE(y,y_hat);
            if error < erbd
                N = N + 1;
%                 fprintf('N %d\n',N);
                out(N,1) = error;
                out(N,2:end) = [b0 bi k dlt p d tau];
            end 

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
