function [alphaK, alphaReturn] = WolfeConditionBisection(f,x0, direction)
%WOLFECONDITIONBISECTION
%   Bisection procedure to select an approximate alpha given a function and
%   a direction over which optimaze on.

bound = 10; % Limit on the selection of the starting alpha (Newton method)
STOP = 0;
c1 = 10e-4;
c2 = 0.9 ; % Working with Newton method

alpha = bound; % Select the initial alpha
alphaK = alpha;
alphaReturn = alpha;
alphaMIN = 0;
alphaMAX = 0;


while STOP == 0
    
    % Fino a quando non rispetto Armijo criterion faccio update di aplha in
    % modo da dimezzare intervallo di ricerca riducendo alphaMAX
    while feval(f, x0 + alpha * direction) > feval(f, x0) + c1 * alpha * grad(f, x0)' * direction 
        alphaMAX = alpha;
        alpha = 0.5*(alphaMIN + alphaMAX);
        alphaReturn = [alphaReturn; alpha];
    end

    % Se alpha trovato al punto soddisfa anche CRVATURE CONDITION allora STOP
    if grad(f, x0 + alpha * direction)' * direction >= c2 * grad(f, x0)' * direction
        alphaK = alpha;
        STOP = 1;
    else
        alphaMIN = alpha;
        
        if alphaMAX == 0
            alpha = 2*alphaMIN;
            alphaReturn = [alphaReturn; alpha];
        else
            alpha = 0.5 * (alphaMIN + alphaMAX);
            alphaReturn = [alphaReturn; alpha];
        end
    end

end
       

end

