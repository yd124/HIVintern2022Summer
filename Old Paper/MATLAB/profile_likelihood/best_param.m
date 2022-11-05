function [x,fval]=best_param(x0,dpi,y) 
    options = optimset('TolFun',1e-6,'TolX',1e-6,'MaxFunEvals', 4000,'MaxIter',4000);
    JJ = @(params) J(params, dpi,y); 
    [x,fval] = fminsearch(JJ,x0,options);
end 
