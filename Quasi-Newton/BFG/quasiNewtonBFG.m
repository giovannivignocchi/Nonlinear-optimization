function [xstar, fstar, counter, error, xks, fks]= quasiNewtonDFP(f, x0, epsilon, maxiterations)

xks = x0';
fks = feavl(f,x0);

counter = 0;
error = 1e300;

% Calcolo H0 utilizzando hessiano SI PUO' FARE MEGLIO ???
H = hes(f,x0);

x = x0;

while error < epsilon || counter < maxiterations
    
    d = - H * grad(x);
    
    % compute an alpha that satisfies weak wolfe condition
    alpha = WolfeConditionBisection(f, x, d);
    
    xk = x + alpha * d;
    
    % Define sk = xk+1 ? xk and yk = ?fk+1 ? ?fk
    sk = xk - x;
    yk = grad(f,xk) - grad(f,x);
    
    H = H + (sk * sk')/(sk'*yk) - (H*(yk*yk')*H)/(yk'*H*yk); 
    
    xks = [xks; xk];
    fks = [fks; feval(f,xk)];
    
    error = norm(grad(f,xk));
    x = xk;
end

xstar = x;
fstar = feval(f,x);

end

