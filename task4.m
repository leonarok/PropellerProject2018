%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Lifting line model with induction factors ------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

load LiftRadius.mat
load('gamma_dist.mat')
wagB=dlmread('wagB.txt'); % J, KT, KQ, eta
geometry= dlmread('Geometry.txt'); % r/R, chord/D, t/D, P/D

itermax=200;
tol=10E-3;

Re=2.4E6;
x0=0.7;
cx0=0.225; % from geometry file
nu=1.19E-6; % from mechanical properties of water
rho=1025;

D=1;
PdivD=1;
AEA0=0.4;
Z=4;
k=100;
damp=0.05;
koeff=0.5;

h=(1-0.167)/(k-1);
x=[0.167:h:1]; % x=r/R
r=x*D/2;
span=D/2*(1-0.167);
phi=atan(PdivD*D./(2*pi*r)); 

chord = interp1(geometry(:,1),geometry(:,2),x,'pchip');
thickness = interp1(geometry(:,1),geometry(:,3),x,'pchip');

J=[0.5 0.6 0.7 0.8 0.9 1.0 1.1];

V= Re*nu/cx0 * 1./sqrt(1+(x0*pi./J).^2);
n= V./(J*D);

C_F = 0.075/(log10(Re)-2)^2; % from ITTC 57'
C_D = 2*C_F*(1+2*thickness./(chord));

C_Lc = interp1(LiftRadius(:,1),LiftRadius(:,2),x,'pchip');


gamma=zeros(length(J),length(x));

U_A=gamma; U_T=gamma; alpha=gamma; Vinf=gamma; C_L=gamma; dTi=gamma;
gamma_new=gamma; gamma_input=gamma;
dQi=gamma; dL=gamma; dT=gamma; dTd=gamma; dKd=gamma; dQ=gamma; dLi=gamma;
dQd=gamma; d_gamma=gamma; beta_i=gamma; i_a=zeros(1,length(x)); i_t=i_a;
twonorm=zeros(1,itermax);

%gamma=2*ones(length(J),length(x));

gamma_0=0.5;
y=-span/2:h*D/2:span/2;
gamma_init=gamma_0*sqrt( 1-( y/(span/2) ).^2 );
for i=1:length(J)
    gamma(i,:)=gamma_init;
end
% for i=1:length(J)
%     gamma(i,:)=interp1(gamma_dist(:,1),gamma_dist(:,i+1),x)/Z;
% end
    
for iter=1:itermax
    for i=1:length(J)
        U_T(i,:)=Z*gamma(i,:)./(2*pi*r);
        for j=1:length(x)
            a=0.5; b=V(i); c=-U_T(i,j).*(2*pi.*r(j)*n(i)-U_T(i,j)/2);
            U_A(i,j)= (-b+sqrt(b^2-4*a*c))/(2*a);
        end
       % beta_i(i,2:end-1)=atan(U_T(i,2:end-1)./U_A(i,2:end-1));
       beta_i(i,1:end)=atan( ( V(i)+0.5*U_A(i,:) )./( 2*pi*r*n(i)-0.5*U_T(i,:) ));
        
       d_gamma(i,:)=gammaderive(D,r,k,gamma(i,:));
        
       for j_fixed=2:length(x)-1
            r0=r(j_fixed);
            
            for j=1:length(x)
                [i_a(j), i_t(j)] = InductionFactors(r(j),...
                    r0,beta_i(i,j),Z);
            end
            U_T(i,j_fixed) = SingularIntegration(r, i_t.*d_gamma(i,:)/(2*pi), r0);
            U_A(i,j_fixed) = SingularIntegration(r, i_a.*d_gamma(i,:)/(2*pi), r0);         
       end
       
        beta_i(i,1:end)=atan(( V(i)+0.5*U_A(i,:) )./(2*pi.*r*n(i)-0.5*U_T(i,:)) );

        alpha(i,2:end-1) = phi(2:end-1)-beta_i(i,2:end-1);
        C_L(i,:) = C_Lc + 2*pi*alpha(i,:);
        
        Vinf(i,:)=sqrt( ( V(i)+0.5*U_A(i,:) ).^2+...
            ( 2*pi*r*n(i)-0.5*U_T(i,:) ).^2);
        
        gamma_new(i,:)=0.5*Vinf(i,:).*chord.*C_L(i,:);
        gamma_input(i,:)=gamma(i,:)+damp*(gamma_new(i,:)-gamma(i,:));
        gamma_input(i,end)=0;gamma_input(i,1)=0;
    end
    twonorm(iter)= (sum(sum(abs(gamma_input-gamma).^2)))^0.5;
    gamma=gamma_input; 
    
    percentdone_print(twonorm(iter)/twonorm(1),tol);
    
    if twonorm(iter)/twonorm(1)<tol
        clc
        fprintf('Convergence criterion met after %d iterations. \n',iter);
        break
    end
end

for i =1:length(J)
dTi(i,:)=rho*gamma(i,:).*(2*pi*r*n(i)-U_T(i,:)/2);
dQi(i,:)=rho*gamma(i,:).*(V(i)+U_A(i,:)/2).*r;
dLi(i,:)=sqrt(dTi(i,:).^2+dQi(i,:).^2/(r.^2));

beta_i(i,:)=atan( (V(i)+0.5*U_A(i,:))/(2*pi.*r*n(i)-0.5*U_T(i,:)) );
dTd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i).*sin(beta_i(i,:));
dQd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i).*r.*cos(beta_i(i,:));

dT(i,:)=dTi(i,:)-dTd(i,:);
dQ(i,:)=dQi(i,:)+dQd(i,:);
end

dT(:,end)=0;
dQ(:,end)=0;

T=Z*trapz(r,dT,2);
Q=Z*trapz(r,dQ,2);

K_T= T'./(rho*n.^2*D^4);
K_Q= Q'./(rho*n.^2*D^5);

eta0=J/(2*pi) .* K_T./K_Q;
figure(1)
hold on
grid minor
title('Thrust, torque & efficiency with induced velocities & induction factors')
xlim([J(1) J(end)])
ylim([0 1.2])
xlabel('J [-]')
ylabel('K_T [-]      10\cdot K_Q [-]      \eta_0 [-]')

plot(J,K_T,'g',wagB(:,1),wagB(:,2),'g--')
plot(J,10*K_Q,'b',wagB(:,1),10*wagB(:,3),'b--')
plot(J,eta0,'r',wagB(:,1),wagB(:,4),'r--')
legend('K_T - lifting line','K_T - Wag. B','10\cdot K_Q - lifting line '...
    ,'10\cdot K_Q - Wag. B','\eta_0 - lifting line','\eta_0 - Wag. B')
