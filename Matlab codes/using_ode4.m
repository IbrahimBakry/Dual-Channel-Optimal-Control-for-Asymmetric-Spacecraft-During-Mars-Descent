% This code: For Aero and mass asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
clc;
clear;
close all;

% Initial Conditions
V0 = 5500; theta0 = -(19*(pi/180)); gamma0= (2*(pi/180)); h0 = 1e5; 
alpha10 = 12*(pi/180); wx10 = 8*(pi/180); thetal0= 0; % thetal (thetaL) NOT theta1 is another equation rather than theta (inclination of irectory)
alpha20 = 12*(pi/180); wx20 = 8*(pi/180); thetat20 = 0;

% Stydy effect of different a & b on wx and alpha
A = [0 1  1   1    ];
B = [1 0.1 0.01 0.005];

Aa = 0.495; Ba = 0.005; Aw = 0.492; Bw = 0.008;

mv = [366, 576, 800];
rv = [1.3, 1.25, 1.6];
Lv = [1.8, 2, 1.8];
Ixv= [135, 443, 506];
Izv= [186, 300, 768];

for i = 1:1 % for A & B 
for j = 2:2 % for KA
epsilon = 0.02; 
% aa = A(i); ba = B(i); aw = A(i); bw = B(i);  % quadratic condition parameters
aa = Aa(i); ba = Ba(i); aw = Aw(i); bw = Bw(i); % Selecting aa,ba,aw,bw values that makes the a(t),wx(t) reach a given interval within a givin time
mzn0 = -0.01; my0f = 0.01; mz0f = 0.01; 
Cx1 = 0.9; Cy1 = 0.9; mz1 = -0.1; 
dzdash = 0.01; dydash = 0.01;
mxw = 0.01; Ixyd = 0.01; Ixzd = 0.01;

y0 = [V0;theta0;gamma0;h0;alpha10;wx10;thetal0;alpha20;wx20;thetat20];
tspan = 0:300; % for discret lines
% tspan = 0:1:300; % for continous lines 

[t,y] = ode45(@(t,y) eqn(t,y,j,aa,ba,aw,bw,epsilon,mzn0,mz1,my0f,...
                         mz0f,mxw,Ixyd,Ixzd),tspan,y0);

V=y(:,1); theta=y(:,2); gamma=y(:,3); h=y(:,4); 
alpha1=y(:,5); wx1=y(:,6); thetal=y(:,7);
alpha2=y(:,8); wx2=y(:,9); thetat2=y(:,10); 

% Maximum deceleration
Cxv = 1.6.*cos(alpha1);
Cya = 1.5.*sin(alpha1);
% a1 = 0.699; a2 = 0.00009; a3 = 47.967; a4 = 0.000426; a5 = a2*a3-a4; a6 = a2*a4;
% Ve = V(1); hs = h;
% a_max = -(Ve^2.*sin(theta(1))./(2*exp(1))).*((a5-a6.*hs)./(a3-a4.*hs));
% n = a_max./g;
rho = (.699.*exp(-0.00009.*h))./(.1921.*((-23.4-0.00222.*h)+273.1));
m = mv(j); r = rv(j); S = pi*r^2; L = Lv(j); Ixd = Ixv(j)/Izv(j);
g0 = 3.72076; Rmars = 3396000;
g = g0.*Rmars.^2./(Rmars+h).^2;
q = 0.5.*rho.*V.^2;
dvdt = Cxv.*(0.5.*rho.*V.^2).*S./m + g.*sin(theta); 
nv = dvdt./g;

% Control laws 
      fw = mxw*q*S^2*L./(Ixd.*V);
     Uwx = -epsilon*(fw+sqrt(fw.^2+aw/bw)).*wx1;
  Ualpha = -epsilon*sqrt(aa/ba).*alpha1;

% Control laws in body coordinate system OXYZ   
Ux = Uwx;
Uy = Ualpha.*0.6; % 0.6 is mid value 
Uz = Ualpha.*cos(gamma);

% Mach Number
gamma_mars = 1.29; % specific heat ratio [dimensionless] 
R_gas_mars = 191.8; % specific gas constant J/kg/K in Metric units
if h >10000, T_mars = -23.4-0.00222.*h +273.1; 
else, T_mars = -32-0.000998.*h +273.1; end
V_sound = sqrt(gamma_mars.*R_gas_mars.*T_mars);

% % Engines ability
% Ix = Ixv(j); Iz = Izv(j); r = rv(j); L = Lv(j); w = sqrt(-mzn0*q*S*L/Iz);
% mux = Ux.*Ix;        % moment about ox
% muy = -Uy.*(2*Iz*w); % moment about oy
% muz = -Uz.*(2*Iz*w); % moment about oz   
% 
% Mux = max(abs(mux)); % Mux = trapz(mux);  % moment [N.m]
% Muy = max(abs(muy)); % Muy = trapz(muy);  % moment [N.m]
% Muz = max(abs(muz)); % Muz = trapz(muz);  % moment [N.m]
% 
% Fux = Mux/(4*r) % N
% Fuy = Muy/(4*r) % N
% Fuz = Muz/(4*r) % N   

%-------------------------------------------------
figure(170); % Hieght with Velocity
    plot(V,h./1000,'LineWidth',3); hold on
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    hold all; grid on; box on; xlabel('Скорость, м/с');ylabel('Высота, км') 
%     ax.YLim = [10000 100000];
    
figure(16); % Deceleration
    plot(nv,h./1000,'LineWidth',3); hold on
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    hold all; grid on; box on; xlabel('Замедление, м/с^2');ylabel('Высота, км') 
%     ax.YLim = [10000 100000];
%     ax.XLim = [1000 t(end)];

figure(1); % alpha
    plot(t,abs(alpha1),'LineWidth',3);  
    ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
    ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
    yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
    yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
    ax.YLim = [0.00 0.25];
    ax.XLim = [0.00 100];
    hold all; grid on; box on; xlabel('t, c');ylabel('Угол атаки \alpha_п, Рад') 

figure(2); % Wx
        plot(t,wx1,'LineWidth',3);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14' '0.16'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
        ax.YLim = [0.00 0.16];
        ax.XLim = [0.00 100];
        hold all; grid on; box on;xlabel('t, c'); ylabel('Угловая скорость \omega_x, 1/c');

figure(2222); % alpha & wx
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
      plot(t,abs(alpha1),'LineWidth',3); hold all
      plot(t,wx1,'LineWidth',3);
      grid on; box on;xlabel('t, c'); ylabel('\omega_x, 1/c, \alpha_п, рад');
      ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
      ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
      ax.YLim = [0.00 0.25];
      ax.XLim = [0.00 100];
ax2 = axes('Position',[0.55 0.55 0.28 0.28]);
      plot(t,abs(alpha1),'LineWidth',3); hold all
      plot(t,wx1,'LineWidth',3);
      grid on; box on;
      ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
      ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
      ax.YLim = [-0.004 0.004];
      ax.XLim = [20 100];
      yticklabels({'-0.004' '-0.002' '0.0' '0.002' '0.004'})
      yticks([-0.004 -0.002 0.0 0.002 0.004])
      xticklabels({'20' '40' '60' '80' '100'})
      xticks([20 40 60 80 100])        
 1       
% figure(22); % ThetaL with control
%         plot(t,thetal,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
% %         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
% %         ax.YLim = [0.00 0.14];
%         hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\thetaL [1/C]');        
        
figure(3); % Ualpha
        plot(t,-abs(Ualpha),'linewidth',3); hold all
        grid on; box on; xlabel('t, c'); ylabel('Функция управления \alpha_п, 1/c')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 50 100 150 200 250 300]) 
%         yticks([-0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.04' '-0.03' '-0.02' '-0.01' '0'})
        %         axis([0 300 -0.1 0])  
        ax.YLim = [-0.03 0]; % ax.YLim = [-0.025 0];
        ax.XLim = [0.00 100];
figure(4); % Uwx
        plot(t,Uwx,'linewidth',3); hold all
        grid on; box on; xlabel('t, c'); ylabel('Функция управления \omega_x 1/C^2');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 25 50 75 100 125 150]) 
%         yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.05' '-0.04' '-0.03' '-0.02' '-0.01' '0'})
%         axis([0 150 -2 0]) 
%         ax.YLim = [-0.06 0]; % ax.YLim = [-0.015 0];
        ax.XLim = [0.00 100];

figure(3333); % Ualpha & Uwx
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
      plot(t,-abs(Ualpha),'LineWidth',3); hold all
      plot(t,Uwx,'LineWidth',3);
      grid on; box on;xlabel('t, c'); ylabel('U\omega_x, 1/c, U\alpha_п, рад');
      ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
      ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
      ax.YLim = [-0.03 0];
      ax.XLim = [0.00 100];
ax2 = axes('Position',[0.55 0.20 0.28 0.28]);
      plot(t,-abs(Ualpha),'LineWidth',3); hold all
      plot(t,Uwx,'LineWidth',3);
      grid on; box on;
      ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
      ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
      ax.YLim = [-0.001 0.001];
      ax.XLim = [20 100];
      yticklabels({'-0.001' '-0.0005' '0.0' '0.0005' '0.001'})
      yticks([-0.001 -0.0005 0.0 0.0005 0.001])
      xticklabels({'20' '40' '60' '80' '100'})
      xticks([20 40 60 80 100])        
1
      
% figure(5); %% Dynamic pressure 
%         plot(t,q,'linewidth',2); ylabel('q') % qv./13
%         grid on; box on; xlabel('t(n) [C]'); hold all
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         xticks([0 50 100 150 200 250 300]) 
%         yticks([0 100 200 300 400 500 600 700])
%         yticklabels({'0', '100', '200', '300', '400', '500', '600', '700'})
%         axis([0 300 0 700]) 

% %% without control
% figure(6); % alpha2 no control
%         plot(t,abs(alpha2),'LineWidth',3);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'-1.0' '-0.75' '-0.50' '-0.25' '0.0' '0.25' '0.50' '0.75' '1.0'})
% %         yticks([-1.0 -0.75 -0.50 -0.25 0.0 0.25 0.50 0.75 1.0]);
% %         ax.YLim = [0.00 0.35];
%         ax.YLim = [0.1 0.3];
%         hold all; grid on; box on;xlabel('t, c');ylabel('Пространственный угол атаки, Рад')  
% figure(7); % wx2 no control
%         plot(t,wx2,'LineWidth',3);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'-0.4' '-0.2' '0.0' '0.2' '0.4' '0.6'})
% %         yticks([-0.4 -0.2 0.0 0.2 0.4 0.6]);
%         ax.YLim = [0.1 0.18];
%         hold all; grid on; box on;xlabel('t, c'); ylabel('Угловая скорость, 1/c');     
% figure(77); % theta (fast phase) no control
%         plot(t,thetat2,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'-0.4' '-0.2' '0.0' '0.2' '0.4' '0.6'})
% %         yticks([-0.4 -0.2 0.0 0.2 0.4 0.6]);
% %         ax.YLim = [0.00 0.14];
%         hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\theta [1/C]');        

% figure(8); % h
%         plot(t,h,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
%         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
%         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('h');  

%---------------------------------------------
% Discretization         
stt = 13;
figure(9); % discret; alpha
        tt = linspace(0,50,stt); aa1 = polyfit(t,alpha1,60); aa1v = polyval(aa1,tt);
        plot(tt,abs(aa1v),'LineWidth',3); hold all;
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
%         yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
        ax.YLim = [0.00 0.25];
        ax.XLim = [0.00 50];
        hold all; grid on; box on;xlabel('t, c');ylabel('Угол атаки \alpha, Рад') 
figure(10); % discret wx 
        tt = linspace(0,50,stt); ww1 = polyfit(t,wx1,20); ww1v = polyval(ww1,tt);
        plot(tt,ww1v,'LineWidth',3); hold all;
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
%         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
        ax.YLim = [0.00 0.15];
        ax.XLim = [0.00 50];
        hold all; grid on; box on;xlabel('t, c');ylabel('Угловая скорость \omega, 1/c');  
stt = 7;
figure(11); % Ualphaz & alphaz
        % to set steps point equal use polyfit
        yyaxis left; hold on;
        tt = linspace(0,50,stt); aa1 = polyfit(t,alpha1,60); aa1v = polyval(aa1,tt);
        plot(tt,abs(aa1v),'LineWidth',3); 
        grid on; box on;xlabel('t, c');ylabel('Угол атаки \alpha, Рад')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.YLim = [0 0.25];
        ax.XLim = [0 50];
        hold off;
        
        yyaxis right; hold on;
        tt = linspace(0,50,stt); UUA = polyfit(t,Ualpha,20); UUV = polyval(UUA,tt);
        stairs(tt,UUV,'linewidth',3); 
        grid on; box on; xlabel('t, c'); ylabel('Функция управления U_\alpha, 1/c')
        ax.XLim = [0 50];
        ax.YLim = [-0.025 0];
        hold off;
        
        ax.XAxis.LineWidth = 4; 
        ax.YAxis(1).LineWidth = 4;
        ax.YAxis(2).LineWidth = 4;
        
figure(12); % Uwxz
        % to set steps point equal use polyfit
        yyaxis right;  hold on;
        tt = linspace(0,50,stt); UUw = polyfit(t,Uwx,20); UUwx = polyval(UUw,tt);
        stairs(tt,UUwx,'linewidth',3);
        grid on; box on; xlabel('t, c'); ylabel('Функция управления U_\omega, 1/c^2');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.YLim = [-0.015 0];
        ax.XLim = [0 50];
        hold off;
        
        yyaxis left; hold on;
        tt = linspace(0,50,stt); ww1 = polyfit(t,wx1,20); ww1v = polyval(ww1,tt);
        plot(tt,ww1v,'LineWidth',3);
        hold all; grid on; box on;xlabel('t, c');ylabel('Угловая скорость \omega, 1/c');  
        ax.YLim = [0 0.15];
        ax.XLim = [0 50];
        hold off;
        
        ax.XAxis.LineWidth = 4; 
        ax.YAxis(1).LineWidth = 4;
        ax.YAxis(2).LineWidth = 4;

stt = 13;
figure(120);    % Ualphaz    
        hold on;
        tt = linspace(0,50,stt); UUA = polyfit(t,Ualpha,20); UUV = polyval(UUA,tt);
        stairs(tt,-abs(UUV),'linewidth',3); 
        grid on; box on; xlabel('t, c'); ylabel('Функция управления U_\alpha, 1/c')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XLim = [0 50];
        ax.YLim = [-0.01 0]; % -0.007 -0.025 -0.07
        hold off;
        ax.XAxis.LineWidth = 4; 
        ax.YAxis.LineWidth = 4;
        
figure(121); % Uwxz
        % to set steps point equal use polyfit
        hold on;
        tt = linspace(0,50,stt); UUw = polyfit(t,Uwx,20); UUwx = polyval(UUw,tt);
        stairs(tt,-abs(UUwx),'linewidth',3);
        grid on; box on; xlabel('t, c'); ylabel('Функция управления U_\omega, 1/c^2');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.YLim = [-0.005 0]; % -0.005 -0.015 -0.05
        ax.XLim = [0 50];
        hold off;
        ax.XAxis.LineWidth = 4; 
        ax.YAxis.LineWidth = 4;
%--------------------------------------------
        
% figure(13); % Gamma
%         plot(t,gamma,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
% %         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
% %         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('Gamma'); 

% figure(14); % Theta "Tang"
%         plot(t,theta,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
% %         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
% %         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('Tang');  

% figure(15); % V
%         plot(t,V,'LineWidth',2);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
% %         yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
% %         ax.YLim = [0.00 100000];
%         hold all; grid on; box on;xlabel('t(n) [C]');ylabel('V'); 

% figure(17); % control functions in OXYZ system, make uw function weight coefficient equal to 1
% subplot(221); % Ux
%         plot(t,-abs(Ux),'linewidth',3); hold all
%         grid on; box on; xlabel('t, с'); ylabel('u_x(t), 1/C^2')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticks([-0.015 -0.010 -0.005 0])
%         yticklabels({'-0.015' '-0.010' '-0.005' '0'})
%         xticks([0 20 40 60 80 100])
%         xticklabels({'0' '20' '40' '60' '80' '100'})
%         ax.YLim = [-0.015 0]; ax.XLim = [0 100];
% subplot(222); % Uy
%         plot(t,-abs(Uy),'linewidth',3); hold all
%         grid on; box on; xlabel('t, c'); ylabel('u_y(t), 1/C');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticks([-0.015 -0.010 -0.005 0])
%         yticklabels({'-0.015' '-0.010' '-0.005' '0'})
%         xticks([0 20 40 60 80 100])
%         xticklabels({'0' '20' '40' '60' '80' '100'})        
%         ax.YLim = [-0.015 0]; ax.XLim = [0 100];
% subplot(223); % Uz
%         plot(t,-abs(Uz),'linewidth',3); hold all
%         grid on; box on; xlabel('t, c'); ylabel('u_z(t), 1/C');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;       
%         yticks([-0.022 -0.015 -0.008 0]); 
%         yticklabels({'-0.022' '-0.015' '-0.008' '0'}) 
%           xticks([0 20 40 60 80 100])
%           xticklabels({'0' '20' '40' '60' '80' '100'})
%         ax.YLim = [-0.022 0]; ax.XLim = [0 100];

% %% Study changing of aa ba aw bw 
% figure(111); % alpha
%     plot(t,abs(alpha1),'LineWidth',3);  
%     ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%     ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %     yticklabels({'0','0.05','0.10','0.15','0.20','0.25'})
% %     yticks([0.0 0.05 0.10 0.15 0.20 0.25]);
% %     ax.YLim = [0.00 0.25];
%     ax.XLim = [0.00 100];
%     hold all; grid on; box on; xlabel('t, c');ylabel('Угол атаки \alpha_п, Рад') 
% figure(222); % Wx
%         plot(t,wx1,'LineWidth',3);
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14' '0.16'})
% %         yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
% %         ax.YLim = [0.00 0.16];
%          ax.XLim = [0.00 100];
%         hold all; grid on; box on;xlabel('t, c'); ylabel('Угловая скорость \omega_x, 1/c');
        
figure(178); % discrete control functions in OXYZ system
subplot(221); % Ux
        % to set steps point equal use polyfit
        tt = linspace(0,100,40); UUx = polyfit(t,Ux,20); UUxv = polyval(UUx,tt);
        stairs(tt,-abs(UUxv),'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([-0.015 -0.01 -0.005 0])
        yticklabels({'-0.015' '-0.01' '-0.005' '0'})
        xticks([0 20 40 60 80 100])
        xticklabels({'0' '20' '40' '60' '80' '100'})        
        ax.YLim = [-0.015 0]; ax.XLim = [0 100];
subplot(222); % Uy
        % to set steps point equal use polyfit
        tt = linspace(0,100,40); UUy = polyfit(t,Uy,20); UUyv = polyval(UUy,tt);
        stairs(tt,-abs(UUyv),'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([-0.015 -0.01 -0.005 0])
        yticklabels({'-0.015' '-0.01' '-0.005' '0'})
        xticks([0 20 40 60 80 100])
        xticklabels({'0' '20' '40' '60' '80' '100'})        
        ax.YLim = [-0.015 0]; ax.XLim = [0 100];       
subplot(223); % Uz
        % to set steps point equal use polyfit
        tt = linspace(0,100,40); UUz = polyfit(t,Uz,20); UUzv = polyval(UUz,tt);
        stairs(tt,-abs(UUzv),'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;       
        yticks([-0.025 -0.02 -0.015 -0.01 -0.005 0])
        yticklabels({'-0.025' '-0.02' '-0.015' '-0.01' '-0.005' '0'})
        xticks([0 20 40 60 80 100])
        xticklabels({'0' '20' '40' '60' '80' '100'})        
        ax.YLim = [-0.025 0]; ax.XLim = [0 100];         

% figure(1780); % discrete control functions in OXYZ system
%         % to set steps point equal use polyfit
%         tt = linspace(0,100,60); UUx = polyfit(t,Ux,20); UUxv = polyval(UUx,tt);
%         stairs(tt,UUxv,'linewidth',3); hold all
%  
%         % to set steps point equal use polyfit
%         tt = linspace(0,100,60); UUy = polyfit(t,Uy,20); UUyv = polyval(UUy,tt);
%         stairs(tt,UUyv,'linewidth',3); hold all
% 
%         % to set steps point equal use polyfit
%         tt = linspace(0,100,60); UUz = polyfit(t,Uz,20); UUzv = polyval(UUz,tt);
%         stairs(tt,UUzv,'linewidth',3); hold all
%         grid on; box on; %xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;       
%         yticks([-0.08 -0.053 -0.026 0])
%         yticklabels({'-0.08' '-0.053' '-0.026' '0'})
%         xticks([0 20 40 60 80 100])
%         xticklabels({'0' '20' '40' '60' '80' '100'})        
%         ax.YLim = [-0.08 0]; ax.XLim = [0 100];         
        
        
% % For Studying Engines Efffect, control functions in OXYZ system 
% figure(190); plot(t,-abs(Ux),'linewidth',3); hold all
% grid on; box on;
% ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
% % ax.YLim = [-0.1 0]; 
% ax.XLim = [0 100];
% %         xticks([0 20 40 60 80 100])
% %         xticklabels({'0' '20' '40' '60' '80' '100'})
% 
% figure(191); plot(t,-abs(Uy),'linewidth',3); hold all
% grid on; box on;
% ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
% % ax.YLim = [-0.1 0]; 
% ax.XLim = [0 100];
% %         xticks([0 20 40 60 80 100])
% %         xticklabels({'0' '20' '40' '60' '80' '100'})
% 
% figure(192); plot(t,-abs(Uz),'linewidth',3); hold all
% grid on; box on;
% ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
% % ax.YLim = [-0.1 0]; 
% ax.XLim = [0 100];
% %         xticks([0 20 40 60 80 100])
% %         xticklabels({'0' '20' '40' '60' '80' '100'})
        
%         st=2;
% figure(20); % For Studying Engines Efffect, control functions in OXYZ system 
%         stairs(nn(1:st:end),Ux(1:st:end),'linewidth',3); hold all
%         stairs(nn(1:st:end),Uy(1:st:end),'linewidth',3); hold all
%         stairs(nn(1:1:end),Uz(1:1:end),'linewidth',3); hold all
%         grid on; box on; 
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
% %         yticks([-0.025 -0.020 -0.015 -0.010 -0.005 0])
% %         yticklabels({'-0.025' '-0.020' '-0.015' '-0.010' '-0.005' '0'})
%         ax.YLim = [-0.08 0];

end
end

function [out] = eqn(t,y,j,aa,ba,aw,bw,epsilon,mzn0,mz1,my0f,mz0f,mxw,Ixyd,Ixzd)
V=y(1); theta=y(2); gamma=y(3); h=y(4); 
alpha1=y(5); wx1=y(6); thetal=y(7);
alpha2=y(8); wx2=y(9); thetat2=y(10);

% «Insight» «Mars Polar Lander» «Mars3»
mv = [366, 576, 800];
rv = [1.3, 1.25, 1.6];
Lv = [1.8, 2, 1.8];
Ixv= [135, 300, 506];
Izv= [186, 443, 768];

g0 = 3.72076; Rmars = 3396000;
r = rv(j); S = pi*r^2; L = Lv(j); m = mv(j); Ix = Ixv(j);Iz = Izv(j);Ixd= Ix/Iz; % Inertia I=Iy=Iz,

% Aerodyamic
Cxv = 1.6.*cos(alpha1);
Cya = 1.5.*sin(alpha1);
 
% Mars Atmospher Density Model
rho = (.699*exp(-0.00009*h))/(.1921*((-23.4-0.00222*h)+273.1));
q = 0.5*(rho)*V^2;
w = sqrt(-mzn0*q*S*L/Iz);
g = g0*Rmars^2/(Rmars+h)^2;

  dVdt = -(Cxv*q*S/m + g*sin(theta));
dthetadt = (- m*g*cos(theta)*(1-V^2/(h+Rmars))/(V*m));
dgammadt = (wx1*cos(alpha1)+tan(theta)*(Cya*q*S*sin(gamma)/(m*V)));
  dhdt = V*sin(theta);     

%---------------------------------------------  
        m1a = -(w^2/mz1)*my0f - Ixzd*wx1^2;
        m2a = -(w^2/mz1)*mz0f + Ixyd*wx1^2;
        mad = sqrt(m1a^2+m2a^2)/w^2;
     theta1 = asin(m1a/mad);

      fw = mxw*q*S^2*L/(Ixd*V);
     Uwx = -epsilon*(fw+sqrt(fw^2+aw/bw))*wx1;
  Ualpha = -epsilon*sqrt(aa/ba)*alpha1;
  
   dwx1dt =  epsilon*mxw*q*S^2*L*wx1/(Ixd*V)+ Uwx;
dalpha1dt = -epsilon*mad*cos(thetal+theta1)*wx1/Ixd + Ualpha;
dthetaldt = (1-0.5*Ixd)*wx1-w;

% NO control
dalpha2dt = -epsilon*mad*cos(thetat2+theta1)*wx2/Ixd;          
   dwx2dt = epsilon*mxw*q*S^2*L*wx2/(Ixd*V);  
dthetat2dt = (1-0.5*Ixd)*wx2-w;
%---------------------------------------------  

out = [dVdt; dthetadt; dgammadt; dhdt; dalpha1dt; dwx1dt; dthetaldt;
       dalpha2dt; dwx2dt; dthetat2dt];
end
