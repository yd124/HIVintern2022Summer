data = readtable('../../data/data1/CHID46.csv');
dpi = data.dpi;
y = data.log_vRNA;
num = size(dpi,1);


h = 0.01;
ti = 0:h:dpi(end); 
init = [10^4 0 10^-3];

mn_mse = Inf;
tic 

N = 0;


params = [0.409*10^-6, 0.233*10^-6, 0.249, 0.775, 14.5*10^3, 0.03, 7];
logV = pred_logV(params, ts, init);

toc 


%% decay function 
function out=b(t,b0,bi,k,tau)
    if t <= tau
        out = b0;
    else
        out = bi+(b0-bi)*exp(-k*(t-tau));
    end
end

function out = pred(params, ts, init)
    
    b0 = params(1);
    bi = params(2);
    k = params(3);
    dlt = params(4);
    p = params(5);
    d = params(6);
    tau = params(7);

    f = @(t,x) [d*(init(1)-x(1))-b(t,b0,bi,k,tau)*x(1)*x(3);...
            b(t,b0,bi,k,tau)*x(1)*x(3)-dlt*x(2);...
            p*x(2)-23*x(3) ];    
    options = odeset('RelTol',1e-4,'AbsTol',1e-6);

    [t,xa] = ode45(f,ts,init,options);
    out = xa;
 
end



%% Mean-squared error
function out = MSE(y,y_hat)
    out = 1/size(y,1) *sum((y-y_hat).^2);
end


% %% Predicting V by viral dynamic model
% 
% function out = pred_logV(params, ts, init)
%     
%     b0 = params(1);
%     bi = params(2);
%     k = params(3);
%     dlt = params(4);
%     p = params(5);
%     d = params(6);
%     tau = params(7);
% 
%     f = @(t,x) [d*(init(1)-x(1))-b(t,b0,bi,k,tau)*x(1)*x(3);...
%             b(t,b0,bi,k,tau)*x(1)*x(3)-dlt*x(2);...
%             p*x(2)-23*x(3) ];    
% 
%     [t,xa] = ode45(f,ts,init);
%     out = log10(xa(:,3));
% 
% end