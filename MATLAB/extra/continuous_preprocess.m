function [dydt, v] = continuous_preprocess(x, y, u, params)
    v = [0, 0];
    dydt = [0, 0, 0, 0, 0, 0];

    if y(1) <= 0
        v(2) = u(2); % / 1 min
    else
        v(2) = max(u(2), y(2));
    end

    dydt(1) = -v(2);
    dydt(2) = v(2);

    if y(6) || u(1) > 0
        v(1) = params.eat_rate;
    end

    dydt(3) = -v(1);

    if y(6) || (~y(6) && v(1) > 0)
        dydt(4) = v(1);
    else
        dydt(4) = 0;
    end

    if v(1) > 0 && ~y(6)
        dydt(6) = true;
        dydt(5) = x(1) + x(2);
    elseif v(1) == 0 && y(6)
        dydt(6) = false;
        dydt(5) = y(5);
    else
        dydt(6) = y(6);
        dydt(5) = y(5);
    end
end