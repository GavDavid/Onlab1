clear;
clf;

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
    Theta(i,1:12) = [1, v_meas(i), s_meas(i),...
    v_meas(i)*s_meas(i), v_meas(i)^2,...
    s_meas(i)^2, s_meas(i)^2*v_meas(i),...
    s_meas(i)*v_meas(i)^2, s_meas(i)^3, v_meas(i)^3, sin(s_meas(i)), cos(s_meas(i))];
end

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
lambda = 1.5;   
n_iter = 10;

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
% zaj és normálás, regularizáció(lambda)



