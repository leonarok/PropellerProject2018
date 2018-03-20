%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%------------- LIFTING LINE METHOD WITHOUT INDUCED VELOCITIES ------------%
%-------------------------------------------------------------------------%
% Calculates torque, thrust and efficiency for a specified propeller      %
% geometry and compares the results with experimental data. This model    %
% does not include the effect of induced velocities                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%-----------------------------------------------------------%
% Load lift coeffiecients, Wagenigen data and geometry data %
%-----------------------------------------------------------%
load LiftRadius.mat;
wagB=dlmread('wagB.txt'); % J, KT, KQ, eta
geometry= dlmread('Geometry.txt'); % r/R, chord/D, t/D, P/D

%------------------------------------%
% Define flow & propeller parameters %
%------------------------------------%
Re=2.4E6;
x0=0.7;
cx0=0.225; % from geometry file
nu=1.19E-6; % from mechanical properties of water
rho=1025;
D=1;
PdivD=1;
AEA0=0.4;
Z=4;

%-----------------------------------------%
% Discretize propeller into foil sections %
%-----------------------------------------%
k=100;
h=(1-0.167)/(k-1);
x=0.167:h:1; % x=r/R
r=x*D/2;
phi=atan(PdivD/pi./x);

%-------------------------------------------------------------%
% Interpolate for chord length and thickness at foil sections %
%-------------------------------------------------------------%
chord = interp1(geometry(:,1),geometry(:,2),x,'pchip');
thickness = interp1(geometry(:,1),geometry(:,3),x,'pchip');

%---------------------------------------------------------------------%
% Initialize vector of advance numbers and find corresponding V and n %
% from eq. in chapter 16 (compendium)                                 %
%---------------------------------------------------------------------%
J=[0.5 0.6 0.7 0.8 0.9 1.0 1.1];
V= Re*nu/cx0 * 1./sqrt(1+(x0*pi./J).^2);
n= V./(J*D);

%-----------------------------------------------------------------%
% Interpolate lift coefficient at zero aoa from the XFoil results %
%-----------------------------------------------------------------%
C_Lc = interp1(LiftRadius(:,1),LiftRadius(:,2),x,'pchip');

%---------------------%
% Initialize matrices %
%---------------------%
beta=zeros(length(J),length(x)); alpha=beta; C_L=beta; Vinf=beta;

%-------------------------------------------------------%
% Calculate thrust and torque on every foil section for %
% every advance number                                  %
%-------------------------------------------------------%
for i =1:length(J)
    
    %--------------------------------------------------------------------%
    % Calculate angle btw velocity and foil and use it to find effective %
    % angle of attack. Find total velocity at every foil section         %
    %--------------------------------------------------------------------%
    beta(i,:) = atan( J(i)./(pi*x) );
    alpha(i,:)=phi-beta(i,:);
    Vinf(i,:)= sqrt( (2*pi*r*n(i)).^2 + V(i).^2 );
    
    %---------------------------------------------------------------%
    % Find lift coeffiecient as the contribution of camber part and %
    % angle of attack part                                          %
    %---------------------------------------------------------------%
    C_L(i,:)= C_Lc(:)' +2*pi*alpha(i,:);
    
    %---------------------------------------------------------------------%
    % Find lift contribution at every foil section, decompose into thrust %
    % and tangential force                                                %
    %---------------------------------------------------------------------%
    dL(i,:)= 0.5*rho*Vinf(i,:).^2.*C_L(i,:)*Z.*chord;
    dTi(i,:)=dL(i,:).*cos(beta(i,:));
    dKi(i,:)=dL(i,:).*sin(beta(i,:));
    
    %---------------------------------------------------------%
    % Find drag coeffient by use of eq. (13.33) in compendium %
    %---------------------------------------------------------%
    Re_c(i,:)=Vinf(i,:).*chord/nu;
    C_F = 0.075./(log10(Re_c(i,:))-2).^2; % from ITTC 57'
    C_D = 2*C_F.*(1+2*thickness./(chord));
    
    %---------------------------------------------------------------------%
    % Find drag contribution at every foil section, decompose into thrust %
    % and tangential force                                                %
    %---------------------------------------------------------------------%
    dD(i,:)=0.5*rho*Vinf(i,:).^2*Z.*chord.*C_D;
    dTd(i,:)=dD(i,:).*sin(beta(i,:));
    dKd(i,:)=dD(i,:).*cos(beta(i,:));
    
    %--------------------------------------------------------------%
    % Find total thrust and torque component at every foil section %
    %--------------------------------------------------------------%
    dT(i,:)=dTi(i,:)-dTd(i,:);
    dQ(i,:)=( dKi(i,:)+dKd(i,:) ).*x*D/2;
end

%--------------------------------------------------%
% Numerical errors at propeller ends are corrected %
%--------------------------------------------------%
dT(:,end)=0;
dQ(:,end)=0;

%--------------------------------------------%
% Integrate with MATLABs trapz() integration %
%--------------------------------------------%
T=trapz(r,dT,2);
Q=trapz(r,dQ,2);

%-------------------------------------------%
% Calculate the thrust and torque coeffient %
%-------------------------------------------%
K_T= T'./(rho*n.^2*D^4);
K_Q= Q'./(rho*n.^2*D^5);

%----------------------%
% Calculate efficiency %
%----------------------%
eta0=J/(2*pi) .* K_T./K_Q;

%--------------------------------------------------------------------%
% Interpolate thrust, torque and effiency from Wag. B for comparison %
%--------------------------------------------------------------------%
K_Twag=interp1(wagB(:,1),wagB(:,2),J);
K_Qwag=interp1(wagB(:,1),wagB(:,3),J);
eta0wag=interp1(wagB(:,1),wagB(:,4),J);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------PLOTTING---------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('rend','painters','pos',[100 100 600 800])
hold on
grid minor
title('Thrust, torque & efficiency without induced velocities')
xlim([J(1) J(end)])
ylim([0 1.25])
xlabel('J [-]')
ylabel('K_T [-]        10\cdot K_Q [-]        \eta_0 [-]')

plot(J,K_T,'g',wagB(:,1),wagB(:,2),'g--')
plot(J,10*K_Q,'b',wagB(:,1),10*wagB(:,3),'b--')
plot(J,eta0,'r',wagB(:,1),wagB(:,4),'r--')
legend('K_T - lifting line','K_T - Wag. B','10\cdot K_Q - lifting line '...
    ,'10\cdot K_Q - Wag. B','\eta_0 - lifting line','\eta_0 - Wag. B')
