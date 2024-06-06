function [dydt, v] = continuous_preprocess(x, y, u, params)
    v = [0, 0];
    dydt = [0, 0, 0, 0, 0, 0];

    if y(1) <= 0
        v(2) = u(2); % / 1 min
    else                
        v(2) = max(u(2), min(y(2), y(1)));  % y1 / dt
    end

    dydt(1) = -v(2);
    dydt(2) = v(2);

    if y(3) > 0 || u(1) > 0
        v(1) = params.eat_rate;
        if y(3) > 0
            v(1) = min(params.eat_rate, y(3));  % y3 / dt
        else
            v(1) = params.eat_rate;
        end
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


function [x_next, y, v] = euler_solve(model, params, x, y_old, u, dt)    
    [dydt, v] = model.preprocess(x,y_old,u,params);
    start_D = y_old(4);

    if u(1) > 0
        start_D = 0;
    end

    y = [u(2) + y_old(1), dydt(2), u(1) + y_old(3), start_D, dydt(5), dydt(6)] + [dydt(1), 0, dydt(3), dydt(4), 0, 0] * dt;
    y = y .* (y >= 0);
    dx_dt = model.step(x,y,v,params);
    x_next = (x + dt * dx_dt);
    x_next = x_next .* (x_next >= 0);
end
