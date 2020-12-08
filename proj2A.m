
clc; clear;
P = 42.66; %kPa
P_bar = 0.4266;

R = 8.314;

%%
x1 = 0:0.05:1;
x2 = 1 - x1;

T_init = 350;
epsilon = 0.0001;

A_met = 11.9673; B_met = 3626.55; C_met = -34.29;
A_water = 11.6834; B_water = 3816.44; C_water = -46.13;




%% Raoult Law
for i = 1:numel(x1)
    y_tot_R(i) = 2;
    T_R(i) = T_init;
    j=1;
    while abs(y_tot_R(i)-1) > epsilon
            
            while j==1
            P1_sat_R(i) = exp(A_met - B_met/(T_R(i) + C_met));
            P2_sat_R(i) = exp(A_water - B_water/(T_R(i) + C_water));

            K1_R(i) = P1_sat_R(i)/P_bar;
            K2_R(i) = P2_sat_R(i)/P_bar;
            
            y_tot_R(i) = K1_R(i) * x1(i) + K2_R(i) * x2(i);
            j = j+1;
            end
            
        
        if y_tot_R(i) > 1
            T_R(i) = 0.9999 * T_R(i);
            
            P1_sat_R(i) = exp(A_met - B_met/(T_R(i) + C_met));
            P2_sat_R(i) = exp(A_water - B_water/(T_R(i) + C_water));

            K1_R(i) = P1_sat_R(i)/P_bar;
            K2_R(i) = P2_sat_R(i)/P_bar;
            
            y_tot_R(i) = K1_R(i) * x1(i) + K2_R(i) * x2(i);
            
        end
        
        if y_tot_R(i) < 1
            T_R(i) = 1.0001 * T_R(i);
            
            P1_sat_R(i) = exp(A_met - B_met/(T_R(i) + C_met));
            P2_sat_R(i) = exp(A_water - B_water/(T_R(i) + C_water));

            K1_R(i) = P1_sat_R(i)/P_bar;
            K2_R(i) = P2_sat_R(i)/P_bar;
            
            y_tot_R(i) = K1_R(i) * x1(i) + K2_R(i) * x2(i);
            
        end
        

    end
    y1_R(i) = K1_R(i) * x1(i);
    
end

plot(x1,T_R,'--b',y1_R,T_R,'--r'); hold on;

%% Wilson

% Using values given for lamdba

%lambda_12 = -1.4636e+03 and lambda_21 = 2.8397e+03;

l_12 = -1.4636e+03; l_21 = 2.8397e+03;
v1 = 41.489; v2 = 18.156;



for k = 1:numel(x1)
    y_tot_W(k) = 2;
    T_W(k) = T_init;
    q=1;
    while abs(y_tot_W(k)-1) > epsilon
            
            while q==1
            L_12(k) = (v2/v1) * exp(-l_12/(R*T_W(k)));
            L_21(k) = (v1/v2) * exp(-l_21/(R*T_W(k)));
            Omega(k) = (L_12(k)/(x1(k)+x2(k)*L_12(k))) - (L_21(k)/(x2(k)+x1(k)*L_21(k)));
            
            gamma_1(k) = exp(-log(x1(k) + x2(k)*L_12(k)) + x2(k)*Omega(k));
            gamma_2(k) = exp(-log(x2(k) + x1(k)*L_21(k)) - x1(k)*Omega(k));
            
            P1_sat_W(k) = exp(A_met - B_met/(T_W(k) + C_met));
            P2_sat_W(k) = exp(A_water - B_water/(T_W(k) + C_water));

            K1_W(k) = P1_sat_W(k)*gamma_1(k)/P_bar;
            K2_W(k) = P2_sat_W(k)*gamma_2(k)/P_bar;
            
            y_tot_W(k) = K1_W(k) * x1(k) + K2_W(k) * x2(k);
            q = q+1;
            end
            
        
        if y_tot_W(k) > 1
            T_W(k) = 0.9999 * T_W(k);
            
            L_12(k) = (v2/v1) * exp(-l_12/(R*T_W(k)));
            L_21(k) = (v1/v2) * exp(-l_21/(R*T_W(k)));
            Omega(k) = (L_12(k)/(x1(k)+x2(k)*L_12(k))) - (L_21(k)/(x2(k)+x1(k)*L_21(k)));
            
            gamma_1(k) = exp(-log(x1(k) + x2(k)*L_12(k)) + x2(k)*Omega(k));
            gamma_2(k) = exp(-log(x2(k) + x1(k)*L_21(k)) - x1(k)*Omega(k));
            
            P1_sat_W(k) = exp(A_met - B_met/(T_W(k) + C_met));
            P2_sat_W(k) = exp(A_water - B_water/(T_W(k) + C_water));

            K1_W(k) = P1_sat_W(k)*gamma_1(k)/P_bar;
            K2_W(k) = P2_sat_W(k)*gamma_2(k)/P_bar;
            
            y_tot_W(k) = K1_W(k) * x1(k) + K2_W(k) * x2(k);
            
        end
        
        if y_tot_W(k) < 1
            T_W(k) = 1.0001 * T_W(k);
            
            L_12(k) = (v2/v1) * exp(-l_12/(R*T_W(k)));
            L_21(k) = (v1/v2) * exp(-l_21/(R*T_W(k)));
            Omega(k) = (L_12(k)/(x1(k)+x2(k)*L_12(k))) - (L_21(k)/(x2(k)+x1(k)*L_21(k)));
            
            gamma_1(k) = exp(-log(x1(k) + x2(k)*L_12(k)) + x2(k)*Omega(k));
            gamma_2(k) = exp(-log(x2(k) + x1(k)*L_21(k)) - x1(k)*Omega(k));
            
            P1_sat_W(k) = exp(A_met - B_met/(T_W(k) + C_met));
            P2_sat_W(k) = exp(A_water - B_water/(T_W(k) + C_water));

            K1_W(k) = P1_sat_W(k)*gamma_1(k)/P_bar;
            K2_W(k) = P2_sat_W(k)*gamma_2(k)/P_bar;
            
            y_tot_W(k) = K1_W(k) * x1(k) + K2_W(k) * x2(k);
            
        end
        

    end
    y1_W(k) = K1_W(k) * x1(k);
    
end

plot(x1,T_W,'-b',y1_W,T_W,'-r'); hold on; grid on;

%% Soujanaya et al (2010)

x1_S = 0:0.1:1;
y1_S = [0 0.41357 0.598483 0.702982 0.77136 0.821163 0.861004 0.89567 0.928488 0.962243 1];
T_S = [350.6 340.4 334.3 330.15 327.1 324.8 322.9 321.3 319.85 318.5 317.15];

plot(x1_S,T_S,'db',y1_S,T_S,'dr');hold on;

xlabel('x_1,y_1');
ylabel('T [K]');

legend('Raoult’s','Wilson','Soujanaya');

