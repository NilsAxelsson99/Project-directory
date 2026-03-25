function x_next = rk4_step(f, t, x, h, params)
    k1 = f(t, x, params);
    k2 = f(t + h/2, x + (h/2) * k1, params);
    k3 = f(t + h/2, x + (h/2) * k2, params);
    k4 = f(t + h,   x + h * k3, params);
    x_next = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end