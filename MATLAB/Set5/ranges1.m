data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);

tic
% out = simulate(10,num,dpi,y,Inf);
% erbd = min(out(:,1));

%0.1*num
out = simulate(1000,dpi,y,0.5);


tb = array2table(out,...
            'VariableNames', ...
            {'J','b0','bi','k','dlt','p','d','tau'});
writetable(tb, 'sim100v2.csv');

toc




%% find the rough minimum error 
function out = simulate(max_N,dpi,y,erbd)
    
    h = 0.01;
    ti = 0:h:dpi(end);

    out = zeros(max_N,8);
    N = 0;
    n = 0;
    
    init = [10^4 0 10^-3];

    while N < max_N
        n = n + 1;
%         fprintf('n %d\n',n);

        %% Unif(a,b) -> a+(b-a)*rand
        %% N(mu,sigma) -> mu+sigma*randn
    
       
        
%         log_bi = -7+2*rand;     
%         log_b0 = -7+2*rand;
%         
% %         log_bi = -6+2*rand;     
% %         log_b0 = -7+1*rand; 
% 
%         if log_bi >= log_b0
%             continue 
%         end
% 
%         log_p = 4*rand;      
%         log_dlt = -1*rand;
%         log_k = 2*rand;   
%         log_d = -1*rand; 
% 
% 
% 
% %         log_k = 2*rand;   
% %         log_p = 4*rand;      
% %         log_dlt = -1*rand;
% %         
% %         log_d = -0.8+0.6*rand; 
% 
%         bi = 10^log_bi;
%         b0 = 10^log_b0;
%         k = 10^log_k;
%         p = 10^log_p;
%         dlt = 10^log_dlt;
%         
%         d = 10^log_d;

%v1

%         b0 = 10^-4*rand;
%         bi = 10^-4*rand;
%         if bi >= b0
%             continue 
%         end 
%         k = 100*rand;
%         p = 5*10^4*rand;
%         dlt = 5*rand;
%         d = 1*rand;

        b0 = 10^-7+(10^-5-10^-7)*rand;
        bi = 10^-9 + (b0 - 10^-9)*rand;
        if bi >= b0
            continue 
        end 
        k = 100*rand;
        p = 2*10^4*rand;
        dlt = 1*rand;
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
            b(t,b0,bi,k,tau)*x(1)*x(3)-10^dlt*x(2);...
            10^p*x(2)-23*x(3) ];    
%     options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [t,xa] = ode45(f,ti,init);
    out = xa;

end




%% Mean-squared error
function out = MSE(y,y_hat)
    out = 1/size(y,1) *sum((y-y_hat).^2);
end
