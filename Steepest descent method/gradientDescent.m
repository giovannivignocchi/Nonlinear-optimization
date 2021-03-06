function [xstar, fstar, counter, error, xks, fks] = gradientDescent(f, x0, epsilon, maxiterations)

  xks = [x0'];
  fks = [feval(f,x0)];

  xk = x0;
  counter = 0;
  error = 1e300;
  
  %Non usiamo dal loop finche non  eccediamo numero max di iterazioni oppure raggiungiamoo un errore ragionevolmete piccolo.
  while error > epsilon && counter < maxiterations

    counter = counter + 1;
    d = - grad(f,xk);
    alpha = fminsearch(@(a) feval(f,xk + a*d), 0.0); %exact 1-d opt
    xk = xk + alpha*d;
	
	% stiamo assumendo di avere una funzione grad che computa il gradiente di f in Xk.
    error = norm(grad(f,xk));
    %aggiorniamo la nostra valutazione di f nel punto Xk cos� trovato.
	fk = feval(f,xk);
    
    xks = [xks; xk'];
    fks = [fks; fk];

  end
  
  xstar = xk;
  fstar = feval(f,xk);

end %of function