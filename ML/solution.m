g = @(t) 9.81 + 0.05 * sin(2 * pi * t);
F = @(t, X) [X(3) * X(4) * g(t) / 25; -X(3)^2 * g(t) / 25; X(1); X(2)];

h = 0.001;
tspan = [0, 100];

X0 = [-4 / 5; -3 / 5; 3; -4];

sol = [];
t = tspan(1):h:tspan(end);

sol(:, 1) = t;
sol(1, 2:5) = X0;
sol(1, 6) = 5;
sol(1, 7) = X0(1) * X0(3) + X0(2) * X0(4);

Xh = X0;
for idx = 2:size(t, 2)
    Xh = rk4(F, t(idx), h, Xh);
    % normalization
    eR = Xh(3:4) / sqrt(Xh(3)^2 + Xh(4)^2);
    Xh(3:4) = 5 * eR;
    % orthogonalization
    eV = Xh(1:2) / sqrt(Xh(1)^2 + Xh(2)^2);
    sphi = eR(1)*eV(1) + eR(2)*eV(2);
    cphi = sqrt(1 - sphi^2);
    tmp = Xh(1) * cphi - sign(Xh(1)) * Xh(2) * sphi;
    Xh(2) = sign(Xh(1)) * Xh(1) * sphi + Xh(2) * cphi;
    Xh(1) = tmp;

    sol(idx, 2:5) = Xh;
    sol(idx, 6) = sqrt(Xh(3)^2 + Xh(4)^2);
    sol(idx, 7) = Xh(1) * Xh(3) + Xh(2) * Xh(4);
end