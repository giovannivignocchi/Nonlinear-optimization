% Simuliamo un' iterazione del Newton method usando weak Wolfe condition
clear;

f = @(x) 3*x(1)^2 + 7*x(2)^2 - 4*x(1) + 9*x(2) - 6;
x1 = -20 + 20*rand;
x2 = -20 + 20*rand;
x0 = [x1 ; x2];

H = hes(f, x0);
G = grad(f, x0);
direction = - inv(H)*G;

[X,Y] = meshgrid(-20:.1:20, -20:.5:20);
Z = zeros(size(X,1),size(X,2));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = feval(f,[X(i,j);Y(i,j)]);
    end
end

figure; contour(X,Y,Z);
hold on; plot(x0(1,1), x0(2,1), 'r*');

[alphak, alphaReturnerd] = WolfeConditionBisection(f, x0, direction);

lim = axis;
axis(lim);

for i=1:size(alphaReturnerd) - 1
    xK = [x0(1,1) + alphaReturnerd(i,1) * direction(1,1); x0(2,1) + alphaReturnerd(i,1) * direction(2,1)];
    hold on; 
    plot(xK(1), xK(2), 'b*');
end

% Plottiamo il punto per il valore di alpha selezionato in colore verde.
hold on; 
plot(x0(1,1) + alphaReturnerd(end,1) * direction(1,1), x0(2,1) + alphaReturnerd(end,1) * direction(2,1), 'g*');
