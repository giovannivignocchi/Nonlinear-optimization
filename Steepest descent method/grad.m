% Gradient of a function at a point
function gradf = grad(f, x)
  h = 0.0001;
  n = length(x);
  gradf = zeros(n,1);
  for i = 1:n
    % Ad ogni iterazione vado a calcolare gradiente nella i-esima direzione
    % cme differenza fra f(x + delta) - f(x) / h, dove delta è uno spostamento
    % infinitesimo nella i-esima direzione e h è il valore di delta.
    delta = zeros(n, 1); delta(i) = h;
    gradf(i) = ( feval(f, x+delta) - feval(f,x) ) / h;
  end
end %of function
