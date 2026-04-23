clear;
clf;
close all;

%% Constants
m = 1;
k = 4;
c = 0.8;
dt = 1/100; % [s]
t = 0:dt:30; % [s]
s_ideal(1) = 1;
v_ideal(1) = 0;

%% Equations
% F = m * a
% -k * x - c * v = m  * a
% -k / m * x - c / m * v = a

for i = 1 : length(t) - 1
    a_ideal(i) = -k / m * s_ideal(i) - c / m * v_ideal(i);
    v_ideal(i+1) = v_ideal(i) + a_ideal(i) * dt;
    s_ideal(i+1) = s_ideal(i) + v_ideal(i) * dt;
end
a_ideal(length(t)) = -k / m * s_ideal(length(t)) - c / m * v_ideal(length(t));

%% SINDy

s_meas = s_ideal;

noise_level = 0.0001;    
gauss = noise_level * randn(size(s_meas));

s_meas = s_meas + gauss;

v_meas = diff(s_meas);
v_meas(length(v_meas) + 1) = v_meas(length(v_meas));
v_meas = v_meas / dt;

a_meas = diff(v_meas);
a_meas(length(a_meas) + 1) = a_meas(length(a_meas));
a_meas = a_meas / dt;


for i = 1 : length(a_ideal)
    Theta(i,1:6) = [1, v_meas(i), s_meas(i),...
    v_meas(i)*s_meas(i), v_meas(i)^2,...
    s_meas(i)^2];
end

Theta_functions = {'1', 'v', 's', 'v*s', 'v^2', 's^2'};

RMS_theta = ones(length(a_ideal), 1) * rms(Theta);
Theta_norm = Theta ./ RMS_theta;

for i = 1 : length(s_meas)
    x(i, 1) = s_meas(i);
    x(i, 2) = v_meas(i);
end

for i = 1 : length(v_meas)
    x_der(i, 1) = v_meas(i);
    x_der(i, 2) = a_meas(i);
end

Xi = inv(Theta_norm'*Theta_norm) * Theta_norm' * x_der(:, 2);

% Sparse regression 
lambda = 0.15;   
n_iter = 100;

for i = 1:n_iter
   
    Xi(abs(Xi) < lambda) = 0;
    
    remains = (Xi ~= 0);
  
    Theta_new = Theta_norm(:, remains);
    Xi_red = inv(Theta_new' * Theta_new) * Theta_new' * x_der(:, 2);
    
    Xi_new = zeros(size(Xi));
    Xi_new(remains) = Xi_red;
    
    Xi = Xi_new;
end

Xi = Xi ./ RMS_theta(1,:)';

% for i = 1:length(Xi)
%     if Xi(i) ~= 0
%         fprintf('%6s tag együtthatója: %8.4f\n', Theta_functions{i}, Xi(i));
%     end
% end

fprintf('\nIdeal equation without noise:\n');
fprintf('a = (%.6f)*v + (%.6f)*s ', (-c / m), (-k / m));
fprintf('\n');

fprintf('\nEquation recreated with  from the noisy measurement:\n');
if Xi(1) ~= 0
    fprintf('a = %.6f', Xi(1));
else
    fprintf('a = ');
end
for i = 2:length(Xi)
    if Xi(i) > 0
        fprintf(' + %.6f*%s', Xi(i), Theta_functions{i});
    elseif Xi(i) < 0
        fprintf(' - %.6f*%s', abs(Xi(i)), Theta_functions{i});
    end
end
fprintf('\n');

% Recreation
s_rec = zeros(1, length(t));
v_rec = zeros(1, length(t));
a_rec = zeros(1, length(t));

s_rec(1) = s_meas(1);
v_rec(1) = v_meas(1);
for i = 1:(length(t) - 1)
    Theta_rec(i,1:6) = [1, v_rec(i), s_rec(i),...
    v_rec(i)*s_rec(i), v_rec(i)^2,...
    s_rec(i)^2];

    a_rec(i) = Theta_rec(i, :) * Xi;
    v_rec(i+1) = v_rec(i) + a_rec(i) * dt;
    s_rec(i+1) = s_rec(i) + v_rec(i) * dt;

end
%% Plot
subplot(3,1,1);
grid on;
grid minor;
hold on;
title 'Distance'
xlabel('[s]');
ylabel('[m]')
plot(t, s_meas, 'r--', 'LineWidth', 0.5);
plot(t, s_ideal, 'b--', 'LineWidth', 0.5);
plot(t, s_rec, 'k', 'LineWidth', 0.5);
legend('Measured', 'Ideal', 'SINDy from measured');

subplot(3,1,2);
grid on;
grid minor;
hold on;
title 'Velocity'
xlabel('[s]');
ylabel('[m/s]')
plot(t, v_meas, 'r--', 'LineWidth', 0.5);
plot(t, v_ideal, 'b--', 'LineWidth', 0.5)
plot(t, v_rec, 'k', 'LineWidth', 0.5);
legend('Measured', 'Ideal', 'SINDy from measured');

subplot(3,1,3);
grid on;
grid minor;
hold on;
title 'Acceleration';
xlabel('[s]');
ylabel('[m/s^2]')
plot(t, a_meas, 'r--', 'LineWidth', 0.5);
plot(t, a_ideal, 'b--', 'LineWidth', 0.5);
plot(t, a_rec, 'k', 'LineWidth', 0.5);
legend('Measured', 'Ideal', 'SINDy from measured');

