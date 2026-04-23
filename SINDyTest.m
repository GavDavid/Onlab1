%% SINDy

dt = 0.1;
t = 1:dt:100;

for i = 1:length(t)
    a(i) = 9.81;
end


for i = 1 : length(a)
    Theta(i,1:10) = [1, v(i), s(i), v(i)*s(i), v(i)^2,...
    s(i)^2, s(i)^2*v(i), s(i)*v(i)^2, s(i)^3, v(i)^3];
end

for i = 1 : length(s)
    x(i, 1) = s(i);
    x(i, 2) = v(i);
end

for i = 1 : length(v)
    x_der(i, 1) = v(i);
    x_der(i, 2) = a(i);
end

Xi= inv(Theta'*Theta) * Theta' * x_der(:, 2);

a_SINDy = Theta * Xi;
v_SINDy = 0;
s_SINDy = 0;
for i = 1 : length(t) - 1
    v_SINDy(i+1) = v_SINDy(i) + a_SINDy(i) * dt;
    s_SINDy(i+1) = s_SINDy(i) + v_SINDy(i) * dt;
end

%% Plot
subplot(3,1,1);
plot(t, a_SINDy, 'b-', t, a, 'r--');
title("Acceleration [m/s^2]")
grid on;
grid minor;
legend("SINDy", "Ideal");
subplot(3,1,2);
plot(t, v_SINDy, 'b-', t, v, 'r--');
title("Velocity [m/s]")
grid on;
grid minor;
legend("SINDy", "Ideal");
subplot(3,1,3);
plot(t, s_SINDy, 'b-', t, s, 'r--');
title("Distance [m]")
grid on;
grid minor;
legend("SINDy", "Ideal");
