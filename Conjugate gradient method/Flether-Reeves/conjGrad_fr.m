% Fletcher-Reeves conjugate gradient method
function [xstar, fstar, counter, error, xks, fks] = conjGrad_fr(f, x0, epsilon, maxiterations)

 % remaining code
  xks = [x0'];
  fks = [feval(f,x0)];
  
  x = x0;
  gradF = grad(f,x);
  d = - gradF;
  
  counter = 0;
  error = 1e300;
  
  %Non usiamo dal loop finche non  eccediamo numero max di iterazioni oppure raggiungiamoo un errore ragionevolmete piccolo.
  while error > epsilon && counter < maxiterations

    counter = counter + 1;
     
    alpha = fminsearch(@(a) feval(f,x + a*d), 0.0);
    
    x = x + alpha * d;
    
    xks = [xks; x'];
    fks = [fks; feval(f,x)];

    error = norm(grad(f,xk));
    
    % Devo memorizzare il balore del gradiente alla precedente iterazione per calcolare Bk
    gradFp = gradF;
    
    % Calcolo gradiente nel punto così trovato
    gradF = grad(f,x);
    
    % Fletcher-Reeves update
    Beta = norm(gradf) / norm(gradFp);
    
    d = - gradf + Beta * d;
    
    xk = xk_new;

  end
  
  xstar = x;
  fstar = feval(f,x);

end %of function

