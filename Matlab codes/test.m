% clc;
% clear;
% close all;

% Atmospher Characteristics
g = 3.9; rho = 0.019;    % Atmospher
r = 1.25; S = pi*r^2; L = 2;             % Sizing
m = 576; Ix = 270; Iz = 443; Ixd= Ix/Iz; % Inertia

% Assumptions
Cxv = 0.04; 
mxw = 0.0001;

% Controller
ba = 1; aw = 0.8; bw = 1;  % these values changes the alphaz, wxz solutions    
hv = 0.5;    % Step of speed equation. 
ht = 0.5;    % Step of theta equation.
hh = 0.5;    % Step of height. 

% Initial Conditions
alpha(1) = 0.32; wx(1) = 0.32;
alphaz(1)= 0.32; wxz(1) = 0.32;
V(1) = 3500; theta(1) = -0.017; h(1) = 1e5;

alpha2(1) = 0.32; wx2(1) = 0.32;
alphaz2(1)= 0.32; wxz2(1) = 0.32;

k = 0;
for i = 0:10:400
          k = k + 1; ii(k) = i;
          q = 0.5*rho*V(k)^2;
% hw = ha condition:
         aa = ba*((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw); 
         Ka = sqrt(aa/ba); 
         ha = 0.3/Ka; % Step of alpha equation.
         ha2(k) = ha;
% Using Euler           
     V(k+1) = V(k) - hv*(Cxv*q*S/m + g*sin(theta(k)));
 theta(k+1) = theta(k) - ht*g*sin(theta(k))/V(k);
     h(k+1) = h(k) + hh*V(k)*sin(theta(k));
         Kw = sqrt((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw); hw = 0.3/Kw; hw2(k)=hw; % Step of angular velocity equation.
 alpha(k+1) = alpha(k)*(1-ha*Ka);
    wx(k+1) = wx(k)*(1-hw*sqrt((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw));
% Using Inverse Z:
alphaz(k+1) = alphaz(k)*(-Ka)^k;
   wxz(k+1) = wxz(k)*(-Kw)^k;
end

k = 0;
for i = 0:10:400
          k = k + 1; ii(k) = i;
          q = 0.5*rho*V(k)^2;
% hw = ha condition:
         aa = 0.8; 
         Ka = sqrt(aa/ba); 
         ha = 0.3/Ka; % Step of alpha equation.
         ha1(k)=ha;
% Using Euler           
     V(k+1) = V(k) - hv*(Cxv*q*S/m + g*sin(theta(k)));
 theta(k+1) = theta(k) - ht*g*sin(theta(k))/V(k);
     h(k+1) = h(k) + hh*V(k)*sin(theta(k));
     
         Kw = sqrt((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw); hw = 0.3/Kw; hw1(k)=hw; % Step of angular velocity equation.
 alpha2(k+1) = alpha2(k)*(1-ha*Ka);
    wx2(k+1) = wx2(k)*(1-hw*sqrt((mxw*q*S^2*L/(Ixd*V(k)))^2+aw/bw));
% Using Inverse Z:
alphaz2(k+1) = alphaz2(k)*(-Ka)^k;
   wxz2(k+1) = wxz2(k)*(-Kw)^k;
end

% Charting
figure; %plot(ii,alpha(1:end-1)); 
        hold on; title('\alpha'); grid on; box on; xlabel('t[s]'); ylabel('\alpha [rad]')
        stairs(ii,abs(alphaz(1:end-1)),'linewidth',2); 
        stairs(ii,abs(alphaz2(1:end-1)),'linewidth',2);  hold off; 
%         legend('Result of Euler Disc.','Result of Z transform')
figure; %plot(ii,alpha(1:end-1)); 
        hold on; title('\alpha'); grid on; box on; xlabel('t[s]'); ylabel('\alpha [rad]') 
        stairs(ii,abs(wxz(1:end-1)),'linewidth',2);
        stairs(ii,abs(wxz2(1:end-1)),'linewidth',2); hold off; 
%         legend('Result of Euler Disc.','Result of Z transform')
