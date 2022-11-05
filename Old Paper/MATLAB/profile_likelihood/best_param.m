function [x,fval]=best_param(x0,dpi,y) 
    options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxFunEvals', 4000,'MaxIter',4000);
    JJ = @(params) J(params, dpi,y); 
    [x,fval] = fminsearch(JJ,x0,options);
end 

%% Least Squared Error
function out=J(params,dpi,y)
    y_hat = pred(params,dpi);
    out =  1/size(y,1) *sum((y-y_hat).^2);
end