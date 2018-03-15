%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Lifting line model with induced velocities -----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

load LiftRadius.mat
wagB=dlmread('wagB.txt'); % J, KT, KQ, eta
geometry= dlmread('Geometry.txt'); % r/R, chord/D, t/D, P/D

itermax=1000;
tol=10E-6;

Re=2.4E6;
x0=0.7;
cx0=0.225; % from geometry file
nu=1.19E-6; % from mechanical properties of water
rho=1025;

D=1;
PdivD=1;
AEA0=0.4;
Z=4;
k=600;
damp=0.05;

h= (1-0.167)/(k-1);
x=0.167:h:1; % x=r/R
r=x*D/2;
phi=atan(PdivD*D./(2*pi*r)); 

chord = interp1(geometry(:,1),geometry(:,2),x,'pchip');
thickness = interp1(geometry(:,1),geometry(:,3),x,'pchip');

J=[0.5 0.6 0.7 0.8 0.9 1.0 1.1];

V= Re*nu/cx0 * 1./sqrt(1+(x0*pi./J).^2);
n= V./(J*D);

C_F = 0.075/(log10(Re)-2)^2; % from ITTC 57'
C_D = 2*C_F*(1+2*thickness./(chord*Z));

C_Lc = interp1(LiftRadius(:,1),LiftRadius(:,2),x,'pchip');

% guess circulation
gamma=zeros(length(J),length(x));
U_A=gamma; U_T=gamma; alpha=gamma; Vinf=gamma; C_L=gamma; dTi=gamma;
dQi=gamma; dL=gamma; dT=gamma; dTd=gamma; dKd=gamma; dQ=gamma; dLi=gamma;
dQd=gamma; gamma_new=gamma; gamma_input=gamma; beta_i=gamma;
twonorm=zeros(1,itermax);

span=D/2*(1-0.167);
gamma_0=10;
y=-span/2:h*D/2:span/2;
gamma_init=gamma_0*sqrt( 1-( y/(span/2) ).^2 );
for i=1:length(J)
    gamma(i,:)=gamma_init;
end

for iter=1:itermax
    for i=1:length(J)
        U_T(i,:)=gamma(i,:)./(2*pi*r);
        for j=1:length(x)
            a=0.5; b=V(i); c=-U_T(i,j).*(2*pi.*x(j)*(D/2)*n(i)-U_T(i,j)/2);
            U_A(i,j)= (-b+sqrt(b^2-4*a*c))/(2*a);
        end
        beta_i(i,2:end-1)=atan(U_T(i,2:end-1)./U_A(i,2:end-1));
        alpha(i,2:end-1) = phi(2:end-1)-beta_i(i,2:end-1);
        C_L(i,:) = C_Lc + 2*pi*alpha(i,:);
        Vinf(i,:)= sqrt( (2*pi*r*n(i)).^2 + V(i).^2 );
        gamma_new(i,:)=0.5*Vinf(i,:)*Z.*chord.*C_L(i,:);
        gamma_input(i,:)=gamma(i,:)+damp*(gamma_new(i,:)-gamma(i,:));
    end
    twonorm(iter)= (sum(sum(abs(gamma_input.^2-gamma.^2))))^0.5;
    gamma=gamma_input;
    if twonorm(iter)/twonorm(1)<tol
        fprintf('Convergence criterion met after %d iterations. \n',iter);
        break
    end
end

beta_i=zeros(length(J),length(x));
for i =1:length(J)
dTi(i,:)=rho*gamma(i,:).*(2*pi*r*n(i)-U_T(i,:)/2);
dQi(i,:)=rho*gamma(i,:).*(V(i)+U_A(i,:)/2).*r;
dLi(i,:)=sqrt(dTi(i,:).^2+dQi(i,:).^2/(r.^2));

beta_i(i,:)=atan( (V(i)+0.5*U_A(i,:))/(2*pi.*r*n(i)-0.5*U_T(i,:)) );
dTd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i)*Z.*sin(beta_i(i,:));
dQd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i)*Z.*r.*cos(beta_i(i,:));

dT(i,:)=dTi(i,:)-dTd(i,:);
dQ(i,:)=dQi(i,:)+dQd(i,:);
end

dT(:,end)=0;
dQ(:,end)=0;

T=trapz(x,dT,2);
Q=trapz(x,dQ,2);

K_T= T'./(rho*n.^2*D^4);
K_Q= Q'./(rho*n.^2*D^5);

eta0=J/(2*pi) .* K_T./K_Q;

% Make gamma-dist file for task 4
gamma_dist(:,1)=x';
gamma_dist(:,2:length(J)+1)=gamma';
save('gamma_dist','gamma_dist')

figure(1)
hold on
grid minor
title('Thrust, torque & efficiency without induced velocities')
xlim([J(1) J(end)])
ylim([0 max(10*K_Q)*1.5])
xlabel('J [-]')
ylabel('K_T [-]      10\cdot K_Q [-]      \eta_0 [-]')

plot(J,K_T,'g',wagB(:,1),wagB(:,2),'g--')
plot(J,10*K_Q,'b',wagB(:,1),10*wagB(:,3),'b--')
plot(J,eta0,'r',wagB(:,1),wagB(:,4),'r--')
legend('K_T - lifting line','K_T - Wag. B','10\cdot K_Q - lifting line '...
    ,'10\cdot K_Q - Wag. B','\eta_0 - lifting line','\eta_0 - Wag. B')
