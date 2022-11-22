%% Least Squared Error
function out=J(params,dpi,y)
    y_hat = pred(params,dpi);
    out =  1/size(y,1) *sum((y-y_hat).^2);
end

%% Predicting V by Using Viral Dynamic Model
function out = pred(params,dpi)
    h = 0.01;
    ti = 0:h:dpi(end);
    init = [1E4 0 1E-3];
    params = 10.^params;

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

%% Decay Function 
function out=b(t,b0,bi,k,tau)
    if t <= tau
        out = b0;
    else
        out = bi+(b0-bi)*exp(-k*(t-tau));
    end
end

