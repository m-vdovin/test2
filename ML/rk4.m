function Xh = rk4(F, t, h, X)

    s1 = F(t, X);
    s2 = F(t + h / 2, X + h * s1 / 2);
    s3 = F(t + h / 2, X + h * s2 / 2);
    s4 = F(t + h, X + h * s3);

    Xh = X + h * (s1 + 2 * s2 + 2 * s3 + s4) / 6;

end