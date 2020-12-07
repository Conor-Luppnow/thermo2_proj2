% Raoult's Law
% y_i*P = x_i * P_i^sat

% Using antoine eqn for P_methanol sat
clc; clear;
P = 42.66; %kPa
P_bar = 0.4266;

% For methanol: A = 11.9673 B = 3626.55 C = -34.29
% For water: A = 11.6834 B = 3816.44 C = -46.13



syms x

A_met = 11.9673; B_met = 3626.55; C_met = -34.29;
A_water = 11.6834; B_water = 3816.44; C_water = -46.13;
x1 = 0:0.05:1;
x2 = 1-x1;

for i=1:numel(x1)
    
    eqn = (x1(i)*exp(A_met - B_met/(x + C_met))/P_bar) + (x2(i)*exp(A_water - B_water/(x + C_water))/P_bar) == 1;
    T(i) = solve(eqn,x);
end


for i=1:numel(x1)
    
    P1_sat(i) = exp(A_met - B_met/(T(i) + C_met))/P_bar;
    P2_sat(i) = exp(A_water - B_water/(T(i) + C_water))/P_bar;
    
    y1(i) = x1(i) * P1_sat(i) / P_bar;
    y2(i) = x2(i) * P2_sat(i) / P_bar;
    
end

plot(x1,T,'ob',y1,T,'or');
grid on;

