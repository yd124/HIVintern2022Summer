%% Least Squared Error
function out=J(params,dpi,y)
    y_hat = pred(params,dpi);
    out =  1/size(y,1) *sum((y-y_hat).^2);
end