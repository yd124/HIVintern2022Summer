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
h=0.01;

tic 
V = zeros(dpi(end)/h+1,1);
for i = 1:num
    
    sti = start:h:dpi(i);
    xa = pred(sti,init,b0,bi,k,dlt,p,d,tau);
    
    idx = round(sti/h+1);
    V(idx) = xa(:,3);
    start = dpi(i);
    init = xa(end,:);

end 
toc 


ti = 0:0.01:dpi(end);
tic
xa2 =pred(ti,[10^4 0 10^-3],b0,bi,k,dlt,p,d,tau);
V2 = xa2(:,3);
toc



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
