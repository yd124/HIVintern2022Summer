data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);

tic
% out = simulate(10,num,dpi,y,Inf);
% erbd = min(out(:,1));


out = simulate(1,num,dpi,y,0.5*num);

params = table2array(readtable('first10.csv'));

x0 = params(1,2:end);



%% find the rough minimum error 
function out = simulate(max_N,num,dpi,y,erbd)
    
    h = 0.01;
    out = zeros(max_N,8);
    N = 0;

    n = 0;
    while N < max_N
        n = n + 1;
        fprintf('n %d\n',n);
        
        params = table2array(readtable('first10.csv'));

        x0 = params(1,2:end);

        b0 = x0(1);
        bi = x0(2);
        k = x0(3);
        dlt = x0(4);
        p = x0(5);
        d = x0(6);
        tau = x0(7);
        
        init = [10^4 0 10^-3];
        start = 0;

        ti = 0:0.01:dpi(end);
        xa2 = pred(ti,[10^4 0 10^-3],b0,bi,k,dlt,p,d,tau);
        V2 = xa2(:,3);
        xa3 = pred(0:h:dpi(1),init,b0,bi,k,dlt,p,d,tau); 
        V3 = xa3(:,3);
        
        V = zeros(dpi(end)/h+1,1);
     
        for i = 1:num
    
            sti = start:h:dpi(i);
            xa = pred(sti,init,b0,bi,k,dlt,p,d,tau);
    
            idx = round(sti/h+1);
            V(idx) = xa(:,3);
            start = dpi(i);
            init = xa(end,:);

            disp(init(3)-V2(dpi(i)/h+1));
            disp(init(3)-V3(dpi(1)/h+1));

        end 




   
        plot(ti,log10(V));
        hold on


        scatter(ti,log10(V2));
        scatter(dpi,y);
        hold off 
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
