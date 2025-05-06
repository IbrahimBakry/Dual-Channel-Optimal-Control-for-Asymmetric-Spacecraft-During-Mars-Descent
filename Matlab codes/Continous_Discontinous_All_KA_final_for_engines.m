% This code: for aero and inertial asymetry
% Theta in "V,h,theta" equations is different from theta in "alpha, omega, theta" equations.
% Notes: more hv, less ht: correct he q curve
%        mzn0 = -0.008 makes the alpha, w curves either closer or far from
%        others.
clc;
clear;
close all;

% «Spirit» «Insight» «Schiaparelli» «Mars Polar Lander» «Mars3»
mv = [174];
rv = [1.15];
Lv = [1.5];
Ixv= [90];
Izv= [80];

for j = 1:length(mv)
% Atmospher Characteristics
g = 3.9; rho = 0.019;  % Mars Atmospher
r = rv(j); S = pi*r^2; L = Lv(j);            % Sizing r = 1.25; L = 2
m = mv(j); Ix = Ixv(j); Iz = Izv(j); Ixd= Ix/Iz; % Inertia I=Iy=Iz, m = 576
Rmars = 3396000;

% Assumptions
 Cxv = 0.04;  epsilon = 0.01; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.01;     mz1 = -0.3;  my0f = 0.3;  mz0f = 0.3; 
dzdash = 0.02; dydash = 0.02;
Ixzd = 0.01;     Ixyd = 0.01; 
 mxw = 0.01;

Ba = [0.01 0.005 0.001 0.0005]; 
Bw = [0.01 0.005 0.001 0.0005]; 

for j = 1:length(Ba)
% Controller
aa = 1; ba = Ba(j); aw = 1; bw = Bw(j);  % quadratic condition parameters
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
for t = 0:4:300 % the step has been changed from 2 to 10 to make the wz & alphaz clearer to the eys
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

mua = -0.2*Ualpha*(2*Iz*w); %wx1(1:end-1)
muw = 0.8*Uwx*Ix;
      
Mua = trapz(mua);  % moment [N.m]
Muw = trapz(muw);  % moment [N.m]

Fua = Mua/(4*rv) % N
Fuw = Muw/(4*rv) % N



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

figure(3); % Ualpha snd Uwx
switch j
    case 1
        subplot(221)
        plot(tt,Ualpha,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        ax.YLim = [-0.04 0];        
        plot(tt,Uwx,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        ax.YLim = [-0.04 0];        
    case 2
        subplot(222)
        plot(tt,Ualpha,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.06 0];
        plot(tt,Uwx,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.06 0];
    case 3
        subplot(223)
        plot(tt,Ualpha,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.1 0];
        plot(tt,Uwx,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.1 0];
    case 4
        subplot(224)
        plot(tt,Ualpha,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Ualpha [Рад]')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.15 0];
        plot(tt,Uwx,'linewidth',3); hold all
        grid on; box on; %xlabel('t(n) [C]'); %ylabel('Uomega [1/C]');
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
                ax.YLim = [-0.15 0];
end
end
end