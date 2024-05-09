function window = set_window(window_size, t, x0, ymin1, u0, t_end)   
    window.t_start = t;
    window.t_end = min(t + window_size, t_end);
    window.size = window.t_end - window.t_start;
    window.x0 = x0;
    window.u0 = u0;
    window.ymin1 = ymin1;
end

