clc;
clear;
close all;

% Atmospher Characteristics
g = 3.9; rho = 0.019;    % Atmospher
r = 1.25; S = pi*r^2; L = 2;             % Sizing
m = 576; Ix = 270; Iz = 443; Ixd= Ix/Iz; II = Iz; % Inertia

% Assumptions
 Cxv = 0.04; epsilon = 0.3; % epsilon 0.02 changes the shape of convergence
mzn0 = -0.0001;  mz1 = 0.1;   my0f = 0.2;    mz0f = 0.2; 
 Cx1 = 0.9;      Cy1 = 0.9; dzdash = 0.05; dydash = 0.05;
Ixzd = 0.01; Ixyd = 0.01;

% Assumptions
Cxv = 0.04; 
mxw = 0.00001;

% Controller
aa = 0.8; ba = 1; aw = 0.8; bw = 1;  % these values changes the alphaz, wxz solutions    
Ka = sqrt(aa/ba); 
hv = 0.5;    % Step of speed equation. 
ht = 0.5;    % Step of theta equation.
hh = 0.5;    % Step of height. 
ha = 0.3;
hw = 0.3;

% Initial Conditions
alphag(1) = 0.32; wxg(1) = 0.24;
alpha1(1) = 0.32; wx1(1) = 0.24;
alphaz1(1)= 0.32; wxz1(1) = 0.24;
V(1) = 3500; theta(1) = -0.017; h(1) = 1e5;
Eps = 0.01; ea = 0.1; ew = 0.1;

k = 0;
for i = 0:10:250
          k = k + 1; ii(k) = i;
% Using Euler           
           q = 0.5*rho*V(k)^2;
      V(k+1) = V(k) - hv*(Cxv*q*S/m + g*sin(theta(k)));
  theta(k+1) = theta(k) - ht*g*sin(theta(k))/V(k);
      h(k+1) = h(k) + hh*V(k)*sin(theta(k));     
      
 alpha1(k+1) = alpha1(k)*(1-ea);
    wx1(k+1) = wx1(k)*(1-ew);
    
% control in the given system
  w = sqrt(-mzn0*q*S*L/II);
m1a = -(w^2/mz1)*my0f - Ixzd*wxg(k)^2;
m2a = -(w^2/mz1)*mz0f - Ixyd*wxg(k)^2;
mad = sqrt(m1a^2+m2a^2);
theta1 = asin(m1a/mad);

         fw = mxw*q*S^2*L*wxz1(k)/(Ixd*V(k));
      uw(k) = -(fw+sqrt(fw^2+aw/bw))*wxz1(k);
      ua(k) = -sqrt(aa/ba)*alphaz1(k);
alphag(k+1) = alphag(k) + ha*(epsilon*mxw*q*S^2*L*wxg(k)/(Ixd*V(k)) + epsilon*uw(k));
   wxg(k+1) = wxg(k) + hw*(-epsilon*mad*cos(theta(k)+theta1)*wxg(k)/Ixd + epsilon*ua(k));

% Using Variable Step
if abs(alpha1(k+1)-alpha1(k)) < Eps; ea = ea;
else, ea = 1.4*ea; end

if abs(wx1(k+1)-wx1(k)) < Eps; ew = ew;
else; ew = 1.3*ew; end

iiz(k) = i;
% Using Inverse Z:
          Kw = sqrt((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw);
alphaz1(k+1) = alphaz1(k)*(-Ka)^k;
   wxz1(k+1) = wxz1(k)*(-Kw)^k;
end

% % Charting
% figure; %plot(ii,alpha(1:end-1)); 
%         hold on; title('\alpha'); grid on; box on; xlabel('t[s]'); ylabel('\alpha [rad]')
%         stairs(ii,alpha1(1:end-1),'linewidth',2); 
%         stairs(ii,abs(alphaz1(1:end-1)),'linewidth',2); hold off; 
%         legend('Result of Euler Disc.','Result of Z transform')
% figure; %plot(ii,wx(1:end-1)); 
%         hold on; title('\omega'); grid on; box on; xlabel('t[s]'); ylabel('\omega [1/s]')
%         stairs(ii,wx1(1:end-1),'linewidth',2);
%         stairs(ii,abs(wxz1(1:end-1)),'linewidth',2); hold off; 
%         legend('Result of Euler Disc.','Result of Z transform')
% 
% figure; bar(ii,abs(alphaz1(1:end-1))-alpha1(1:end-1)); title('\alpha Solution accuracy') 
%         grid on; box on; xlabel('t[s]'); ylabel('\delta\alpha [rad]')
% figure; bar(ii,abs(wxz1(1:end-1))-wx1(1:end-1)); title('\omega Solution accuracy') 
%         grid on; box on; xlabel('t[s]'); ylabel('\delta\omega [rad/cek]')
%         

% Charting
figure; %plot(ii,alpha(1:end-1)); 
        hold on; grid on; box on;
        xlabel('t(n) [C]'); 
        ylabel('\alpha [Рад]')
        stairs(ii,alpha1(1:end-1),'g','linewidth',6); 
        stairs(iiz,abs(alphaz1(1:end-1)),'c','linewidth',6); hold off; 
%         stairs(ii,abs(alphag(1:end-1)),'c','linewidth',6); hold off;         
        legend('Result of Euler Disc.....................',...
               'Result of Z transform......,,,,,,,,,,....')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35])
        yticklabels({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30', '0.35'})
        ax.YLim = [0.00 0.35];
        
figure; %plot(ii,wx(1:end-1)); 
        hold on; grid on; box on; 
        xlabel('t(n) [C]'); 
        ylabel('\omega [1/C]');
        stairs(ii,wx1(1:end-1),'g','linewidth',6);
        stairs(iiz,abs(wxz1(1:end-1)),'c','linewidth',6); hold off; 
%         stairs(ii,abs(wxg(1:end-1)),'c','linewidth',6); hold off;         
        legend('Result of Euler Disc......................',...
               'Result of Z transform.....................')
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
        yticks([0.00 0.05 0.10 0.15 0.20 0.25])
        ax.YLim = [0.00 0.25];

%% In the given system        
% figure; 
%         grid on; box on; xlabel('t(n) [C]'); ylabel('\alpha [Рад]')
%         stairs(ii,abs(alphag(1:end-1)),'c','linewidth',6); hold off;         
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;
%         
% figure;  
%         grid on; box on; xlabel('t(n) [C]'); ylabel('\omega [1/C]');
%         stairs(ii,abs(wxg(1:end-1)),'c','linewidth',6);       
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4;

%% The controls
figure; 
        stairs(ii,-abs(ua),'g','linewidth',6); hold off;         
        ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 20; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
        ax.XAxis.LineWidth = 4; ax.YAxis.LineWidth = 4; hold on 
        stairs(ii,-abs(uw),'b','linewidth',6); hold off 
        grid on; box on; xlabel('t(n) [C]'); ylabel('\alpha [Рад]')
        yticks([-0.3 -0.25 -0.20 -0.15 -0.10 -0.05 0])
        yticklabels({'-0.30', '-0.25', '-0.20', '-0.15', '-0.10', '-0.05', '0'})
        legend('U\alpha',...
               'U\omegaz')
%%
% figure; bar(ii,abs(alphaz1(1:end-1))-alpha1(1:end-1));
%         grid on; box on; 
%         xlabel('t(n) [C]'); 
%         ylabel('\delta\alpha [рад]')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
%         yticks([-0.04 -0.03 -0.02 -0.01 0.01])
%         ax.YLim = [-0.04 0.01];

% figure; bar(ii,abs(wxz1(1:end-1))-wx1(1:end-1));
%         grid on; box on; 
%         xlabel('t(n) [C]'); 
%         ylabel('\delta\omega [1/C]')
%         ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
%         ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
%         yticks([-0.01 -0.005 0.005 0.01])
%         ax.YLim = [-0.01 0.01];
        
                
        
        
        