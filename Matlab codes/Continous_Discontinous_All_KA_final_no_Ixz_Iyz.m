% This code: for aero and inertial asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
% Notes: more hv, less ht: correct he q curve
%        mzn0 = -0.008 makes the alpha, w curves either closer or far from
%        others.
clc;
clear;
close all;

% «Spirit» «Insight» «Schiaparelli» «Mars Polar Lander» «Mars3»
mv = [174, 366, 577, 576, 800];
rv = [1.15, 1.3, 1.2, 1.25, 1.6];
Lv = [1.5, 1.8, 1.8, 2, 1.8];
Ixv= [90, 186, 250, 270, 768];
Izv= [80, 135, 195, 443, 506];

for j = 1:length(mv)
% Atmospher Characteristics
g = 3.9; rho = 0.019;  % Mars Atmospher
r = rv(j); S = pi*r^2; L = Lv(j);            % Sizing r = 1.25; L = 2
m = mv(j); Ix = Ixv(j); Iz = Izv(j); Ixd= Ix/Iz; % Inertia I=Iy=Iz, m = 576
Rmars = 3396000;

% Assumptions
 Cxv = 0.04; epsilon = 0.01; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.01;  mz1 = -0.3;  my0f = 0.3;  mz0f = 0.3; 
dzdash = 0.02; dydash = 0.02;
Ixzd = 0;     Ixyd = 0; 
 mxw = 0.01;

% Controller
aa = 1; ba = 0.005; aw = 1; bw = 0.005;  % quadratic condition parameters
hv = 10; ht = 0.001; hh = 12; 
ha = 2;
hw = 2;
% high hv, small ht: forms the right q curve 

% Initial Conditions
alpha1(1) = 0.32; wx1(1) = 0.14; thetal(1) = 0; % thetal (thetaL) NOT theta1 is another equation rather than theta (inclination of irectory)
alpha2(1) = 0.32; wx2(1) = 0.14;
alphaz(1) = 0.32; wxz(1) = 0.14;
V(1) = 3500; theta(1) = -0.017; h(1) = 1e5;

k = 0;
for t = 0:8:300 % the step has been changed from 2 to 10 to make the wz & alphaz clearer to the eys
    k = k + 1;  kk(k) = k; tt(k) = t; 
             
     [rho,~] = marsatmoshper(h(k)); % Mars Atmospher Density Model
           q = 0.5*(rho)*V(k)^2;     RHO(k) = rho; qv(k)=q;
           w = sqrt(-mzn0*q*S*L/Iz);
       
      V(k+1) = V(k) - hv*(Cxv*q*S/m + g*sin(theta(k)));
  theta(k+1) = theta(k) - ht*(cos(theta(k))*(g-V(k)^2/(h(k)+Rmars))/V(k));
      h(k+1) = h(k) + hh*V(k)*sin(theta(k));     

% % Using Euler - control in the given system
        m1a = -(w^2/mz1)*my0f - Ixzd*wx1(k)^2;
        m2a = -(w^2/mz1)*mz0f + Ixyd*wx1(k)^2;
        mad = sqrt(m1a^2+m2a^2)/w^2;
     theta1 = asin(m1a/mad);

         fw = mxw*q*S^2*L/(Ixd*V(k));
     Uwx(k) = -epsilon*(fw+sqrt(fw^2+aw/bw))*wx1(k);
  Ualpha(k) = -epsilon*sqrt(aa/ba)*alpha1(k);
   wx1(k+1) = wx1(k) + hw*(epsilon*mxw*q*S^2*L*wx1(k)/(Ixd*V(k)) + Uwx(k));
alpha1(k+1) = alpha1(k) + ha*(-epsilon*mad*cos(thetal(k)+theta1)*wx1(k)/Ixd + Ualpha(k));
thetal(k+1) = thetal(k) - 0.1*epsilon*(1-0.5*Ixd)*wx1(k)-w;

% NO control
alpha2(k+1) = alpha2(k) + ha*(-epsilon*mad*cos(thetal(k)+theta1)*wx2(k)/Ixd);          
   wx2(k+1) = wx2(k) + hw*(epsilon*mxw*q*S^2*L*wx2(k)/(Ixd*V(k)));  
   
end


% Special for loop for the wxz & alphaz, to make there data longer 
Vz(1) = 3500; thetaz(1) = -0.017; hz(1) = 1e5;
k = 0;
for t = 0:1:600 
    k = k + 1;  kkz(k) = k; ttz(k) = t; 
     [rho,~] = marsatmoshper(hz(k)); % Mars Atmospher Density Model
           q = 0.5*(rho)*Vz(k)^2; 
       
      Vz(k+1) = Vz(k) - hv*(Cxv*q*S/m + g*sin(thetaz(k)));
  thetaz(k+1) = thetaz(k) - ht*(cos(thetaz(k))*(g-Vz(k)^2/(hz(k)+Rmars))/Vz(k));
      hz(k+1) = hz(k) + hh*Vz(k)*sin(thetaz(k));
      
      % Discontinous Using Inverse Z:
          Ka = sqrt(aa/1.05); % ba & bw smaller to move the curve to right, but not equal or smaller then 1
          Kw = sqrt((mxw*q*S^2*L/(Ixd*Vz(k)))^2+aw/1.05);
 alphaz(k+1) = alphaz(k)*(-Ka)^k;
    wxz(k+1) = wxz(k)*(-Kw)^k;
end

% figure;hold all
% subplot(221); plot(tt,h(1:end-1)); hold all
% subplot(222); plot(tt,V(1:end-1)); hold all
% subplot(223); plot(tt,theta(1:end-1)); hold all

% Charting
figure(1); % alpha
        plot(tt,(alpha1(1:end-1)),'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
        yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
        ax.YLim = [0.00 0.35];
        ax.XLim = [0.00 tt(end)];
        hold all; grid on; box on; xlabel('t(n) [C]');ylabel('\alpha [Рад]') 
        
figure(2); % Wx
        plot(tt,wx1(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
        ax.YLim = [0.00 0.14];
        hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\omega [1/C]');

figure(3); % Ualpha
        plot(tt,Ualpha,'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 50 100 150 200 250 300]) 
%         yticks([-0.04 -0.03 -0.02 -0.01 0])
%         yticklabels({'-0.04' '-0.03' '-0.02' '-0.01' '0'})
%         %         axis([0 300 -0.1 0])  
%         ax.YLim = [-0.04 0];
figure(4); % Uwx
        plot(tt,Uwx,'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 25 50 75 100 125 150]) 
        yticks([-0.020 -0.015 -0.010 -0.005 0])
        yticklabels({'-0.020' '-0.015' '-0.010' '-0.005' '0'})
%         axis([0 150 -2 0]) 
        ax.YLim = [-0.02 0];
figure(5); %% Dynamic pressure 
        plot(tt,qv./13,'linewidth',2); ylabel('q')
        grid on; box on; xlabel('t(n) [C]'); hold all
        ax = gca; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4; ax.GridLineStyle = ':';
% %         xticks([0 50 100 150 200 250 300]) 
%         yticks([0 100 200 300 400 500 600 700])
%         yticklabels({'0', '100', '200', '300', '400', '500', '600', '700'})
%         axis([0 300 0 700]) 

% % without control
figure(6); % alpha2 no control
        plot(tt,(alpha2(1:end-1)),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'-1.0' '-0.75' '-0.50' '-0.25' '0.0' '0.25' '0.50' '0.75' '1.0'})
%         yticks([-1.0 -0.75 -0.50 -0.25 0.0 0.25 0.50 0.75 1.0]);
% %         ax.YLim = [0.00 0.35];
% %         ax.XLim = [0.00 tt(end)];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\alpha [Рад]')     
figure(7); % wx2 no control
        plot(tt,wx2(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         yticklabels({'-0.4' '-0.2' '0.0' '0.2' '0.4' '0.6'})
%         yticks([-0.4 -0.2 0.0 0.2 0.4 0.6]);
% %         ax.YLim = [0.00 0.14];
        hold all; grid on; box on;xlabel('t(n) [C]'); ylabel('\omegax [1/C]');
figure(8); % h
        plot(tt,h(1:end-1),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'})
        yticks([0.0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000])
        ax.YLim = [0.00 100000];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('h');  
        
figure(9); % discret alpha
        plot(tt(1:1:end),alpha1(1:1:end-1),'LineWidth',2);hold all
        plot(1:10:252,abs(alphaz(1:2:51)),'LineWidth',2);
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'})
        yticks([0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.345]);
        ax.YLim = [0.00 0.35];
        ax.XLim = [0.00 200];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\alpha [Рад]') 
figure(10); % discret wx euler and z
        plot(tt,wx1(1:end-1),'LineWidth',2); hold all
        plot(1:10:252,abs(wxz(1:2:51)),'LineWidth',2); 
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticklabels({'0.0' '0.02' '0.04' '0.06' '0.08' '0.10' '0.12' '0.14'})
        yticks([0.0 0.02 0.04 0.06 0.08 0.10 0.12 0.14])
        ax.YLim = [0.00 0.14];
        ax.XLim = [0.00 200];
        hold all; grid on; box on;xlabel('t(n) [C]');ylabel('\omega [1/C]');  
figure(11); % Ualphaz
        stairs(tt,Ualpha,'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 50 100 150 200 250 300]) 
        yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0])
        yticklabels({'-0.05' '-0.04' '-0.03' '-0.02' '-0.01' '0'})
        %         axis([0 300 -0.1 0])  
        ax.YLim = [-0.05 0];
        ax.XLim = [0.00 200];
figure(12); % Uwxz
        stairs(tt,Uwx,'linewidth',2); hold all
        grid on; box on; xlabel('t(n) [C]'); ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         xticks([0 25 50 75 100 125 150]) 
        yticks([-0.020 -0.015 -0.010 -0.005 0])
        yticklabels({'-0.020' '-0.015' '-0.010' '-0.005' '0'})
%         axis([0 150 -2 0]) 
        ax.YLim = [-0.02 0];
        ax.XLim = [0.00 200];
 
end