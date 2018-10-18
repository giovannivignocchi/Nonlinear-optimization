function [xstar, vstar] = ipm(c, A, b, mu, x_init, gamma, epsilon)

  OPTIONS = [];
  [m, n] = size(A);
  
  %bring to standard form
  % interior point method for LPs of the form min c'*x : A*x <= b
  
  s = b - A*x_init; % Assegno valore alle slack in funzione del punto iniziale x_init
  x = [ x_init ; s ]; % Aggiungo le slack variables
  
  A = [A eye(m,m)]; % Affianco matrce identità a coefficienti constraints per tenere conto slack variables
  c = [c zeros(1,m)]; % Inserisco coefficiente nullo in obj per le slack variables



  %preparo le matrici che userò nel ciclo while
  grad = zeros(m+n, 1);
  H = zeros(m+n, m+n);
  d = zeros(m+n, 1);
  u = zeros(m, 1);
  
  xks = [x_init'];

  %% iterative method
  while n * mu >= epsilon
    %% compute gradient and Hessian
	for i=1:(n+m)
        grad(i) = c(i) - mu/x(i);
        H(i,i) = mu/x(i)^2;
    end
    %% compute direction d with adapted Newton update
    
    % Compose the KKT Matrix
    N = [ H, A';
          A, zeros(m,m)];
    
    bN = [-grad; zeros(m,1)];
    
    xN = N \ bN;
    
    d = xN(1 : n+m);
    u = xN(n+m+1 : end);
    

    alpha = fminbnd(@(alpha) c*(x + alpha*d) - mu*sum(log(x + alpha*d),1), 0, 1, OPTIONS);
    xstar = x + alpha * d;
    mu = gamma * mu;
    x = xstar;
    xks = [xks; x(1:n)'];
  end
  

%polyhedron_print(A,b); hold on;
plot(xks(:,1), xks(:,2), 'r.');

end % end of function

