g = @(t) 9.81 + 0.05 * sin(2 * pi * t);
F = @(t, X) [-sin(X(2)) * g(t) / 5; X(1)];

h = 0.001;
tspan = [0, 100];

X0 = [-1 / 5; atan2(3, 4)];

ref = [];
t = tspan(1):h:tspan(end);

ref(:, 1) = t;
ref(1, 2:3) = X0;
ref(1, 4) = 5 * sin(X0(2));
ref(1, 5) = -5 * cos(X0(2));

Xh = X0;
for idx = 2:size(t, 2)
    Xh = rk4(F, t(idx), h, Xh);

    ref(idx, 2:3) = Xh;
    ref(idx, 4) = 5 * sin(Xh(2));
    ref(idx, 5) = -5 * cos(Xh(2));
end