function [xstar, fstar, counter, error, flocals, fstars, nstars] = ...
    multiStart(f, n, eps, localeps, maxit, maxlocalit, myLocalOptimAlg)

fstars = []; % estimate of the optimal value of the global objective function
nstars = []; % iteration in which is updated the estimate of the optimal value
flocals = []; % optimal value of the local objective function

bound = 5;

% n = (numberOfVariables - 1) * component(3)
x = rnd(n, bound);
xstar = x;
counter = 0;
termination = 0;

while termination == 0
    
    % Con questa fstar stiamo andando a calcolare la somma degli
    % Square-error che otteniamo con le posizioni degli atomi calcolati
    % all'iterazione precedente.
    %
    % Se l'errore commesso è più piccolo dell'epsilon fornito in input alla
    % funzione multiStart fermiamo l'esecuzione.
    
    fstar = feval(f, xstar);
    if fstar < eps || counter >= maxit
        termination = 1;
        
    else
        % Nuova iterazione aumentiamo il contatore.
        counter = counter + 1;
        
        [xlocal, flocal, count, err] = myLocalOptimAlg(f,x,localeps,maxlocalit);
        flocals = [flocals; flocal];
            
            if flocal < fstar
                xstar = xlocal;
                fstar = flocal;
                error = err;
                %
                fstars = [fstars;fstar];
                nstars = [nstars;counter];
            end
        % Preparo il punto da cui partire nella prossima iterazione.
        x = rnd(n, bound);
    end
end
end