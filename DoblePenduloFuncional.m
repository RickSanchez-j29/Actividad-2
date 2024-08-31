function double_pendulo
    % Parámetros del sistema
    m1 = 0.5; % kg
    m2 = 0.5; % kg
    l1 = 0.3; % m
    l2 = 0.3; % m
    g = 9.81; % m/s^2

    % Condiciones iniciales
    theta1 = 0; % rad
    theta2 = pi/6; % rad
    dtheta1 = 0; % rad/s
    dtheta2 = 0; % rad/s
    y0 = [theta1; dtheta1; theta2; dtheta2];

    % Intervalo de tiempo
    tspan = [0 10]; % segundos

    % Resolver el sistema de ecuaciones diferenciales
    [t, y] = ode45(@(t, y) pendulo_doble(t, y, m1, m2, l1, l2, g), tspan, y0);

    % Graficar los resultados
    figure;
    plot(t, y(:, 1), t, y(:, 3));
    xlabel('Tiempo (s)');
    ylabel('Ángulo (rad)');
    legend('\theta_1', '\theta_2');
    title('Movimiento del doble péndulo');
end

function dydt = pendulo_doble(~, y, m1, m2, l1, l2, g)
    theta1 = y(1);
    dtheta1 = y(2);
    theta2 = y(3);
    dtheta2 = y(4);

    % Ecuaciones del movimiento
    
    delta = theta2 - theta1;
    den1 = (m1 + m2) * l1 - m2 * l1 * cos(delta)^2;
    den2 = (l2 / l1) * den1;

    
    d2theta1 = (m2 * l1 * dtheta1^2 * sin(delta) * cos(delta) + m2 * g * sin(theta2) * cos(delta) + m2 * l2 * dtheta2^2 * sin(delta) - (m1 + m2) * g * sin(theta1)) / den1;
    
    d2theta2 = (-m2 * l2 * dtheta2^2 * sin(delta) * cos(delta) + (m1 + m2) * (g * sin(theta1) * cos(delta) - l1 * dtheta1^2 * sin(delta) - g * sin(theta2))) / den2;

    dydt = [dtheta1; d2theta1; dtheta2; d2theta2];
end
