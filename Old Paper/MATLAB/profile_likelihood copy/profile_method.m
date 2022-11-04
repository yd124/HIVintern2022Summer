


%% Profile Methods
function out = profile_method(xx,sampled_params,num_params, num_samples,dpi,y)
    cd("~/Documents/MATLAB")
    options = optimoptions(@lsqnonlin,'MaxIterations', 500);
    LSE = zeros(num_params,num_samples);
    for ip = 1:num_params
        for is = 1:num_samples
            fp = sampled_params(ip,is);
            ic = [xx(1:ip-1) xx(ip+1:end)];
            fun = @(r) pred([r(1:ip-1) fp r(ip:end)],dpi) - y;
            [x,resnorm] = lsqnonlin(fun,ic,zeros(size(ic)),[],options);
            LSE(ip,is) = resnorm/size(y,1);
        end
    end
    out = LSE;
end 

%% Decay Function 
function out=b(t,b0,bi,k,tau)
    if t <= tau
        out = b0;
    else
        out = bi+(b0-bi)*exp(-k*(t-tau));
    end
end

%% Predicting V by Using Viral Dynamic Model
function out = pred(params,dpi)
    h = 0.01;
    ti = 0:h:dpi(end);
    init = [10^4 0 10^-3];
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

    [t,xa] = ode45(f,ti,init);
    V = xa(:,3);
    try
        out = log10(V(dpi/h+1));
    catch 
        out = zeros(size(dpi));
    end
end

%% Least Squared Error
function out=J(params,dpi,y)
    y_hat = pred(params,dpi);
    out =  1/size(y,1) *sum((y-y_hat).^2);
end
