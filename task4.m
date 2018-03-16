%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-LIFTING LINE METHOD WITH INDUCTION FACTOR METHOD FOR INDUCED VELOCITIES-%
%-------------------------------------------------------------------------%
% Calculates torque, thrust and efficiency for a specified propeller      %
% geometry and compares the results with experimental data. This model    %
% includes the effect of induced velocities through the induction         %
% factor method.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%----------------------------------------------------------------%
% Load lift coeffiecients, circulation distribution from task 4, %
% Wagenigen data and geometry data                               %
%----------------------------------------------------------------%
load LiftRadius.mat;
load('gamma_dist.mat');
wagB=dlmread('wagB.txt'); % J, KT, KQ, eta
geometry= dlmread('Geometry.txt'); % r/R, chord/D, t/D, P/D

%------------------------------%
% Specify iteration parameters %
%------------------------------%
itermax=1000;
tol=10E-4;
damp=0.05;

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
span=D/2*(1-0.167);

%-----------------------------------------%
% Discretize propeller into foil sections %
%-----------------------------------------%
k=60;
h=(1-0.167)/(k-1);
x=0.167:h:1; % x=r/R
r=x*D/2;
phi=atan(PdivD*D./(2*pi*r)); 

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

%---------------------------------------------------------%
% Find drag coeffient by use of eq. (13.33) in compendium %
%---------------------------------------------------------%
C_F = 0.075/(log10(Re)-2)^2; % from ITTC 57'
C_D = 2*C_F*(1+2*thickness./(chord));

%-----------------------------------------------------------------%
% Interpolate lift coefficient at zero aoa from the XFoil results %
%-----------------------------------------------------------------%
C_Lc = interp1(LiftRadius(:,1),LiftRadius(:,2),x,'pchip');

%-------------------------------%
% Initialize matrices & vectors %
%-------------------------------%
gamma=zeros(length(J),length(x));
U_A=gamma; U_T=gamma; alpha=gamma; Vinf=gamma; C_L=gamma; dTi=gamma;
gamma_new=gamma; gamma_input=gamma;
dQi=gamma; dL=gamma; dT=gamma; dTd=gamma; dKd=gamma; dQ=gamma; dLi=gamma;
dQd=gamma; d_gamma=gamma; beta_i=gamma; i_a=zeros(1,length(x)); i_t=i_a;
twonorm=zeros(1,itermax);

%----------------------------------------------------------------%
% Interpolate circulation from the circulation dist. from task 3 %
%----------------------------------------------------------------%
for i=1:length(J)
    gamma(i,:)=interp1(gamma_dist(:,1),gamma_dist(:,i+1),x)/Z;
end

%-------------------------------------------------------------------------%
% Iterate the procedure until convergence is satisfied or maximum allowed %
% number of iterations is reached                                         %
%-------------------------------------------------------------------------%
for iter=1:itermax
    %---------------------------------------------------------%
    % Calculate the new circulation at every foil section for %
    % every advance number                                    %
    %---------------------------------------------------------%
    for i=1:length(J)
        
        %--------------------------------------------------------------%
        % Calculate mean tangential and axial induced velocities using %
        % complete momentum theory at every foil section.              %
        %--------------------------------------------------------------%
        U_T(i,:)=Z*gamma(i,:)./(2*pi*r); % eq. (13.43)
        for j=1:length(x)
            a=0.5; b=V(i); c=-U_T(i,j).*(2*pi.*r(j)*n(i)-U_T(i,j)/2);
            U_A(i,j)= (-b+sqrt(b^2-4*a*c))/(2*a); % eq. (13.28)
        end

        %-----------------------------------------------------------%
        % Calculate tentative hydrodynamic pitch angle and find the %
        % derivative of the circ. dist. over the propeller blade    %
        %-----------------------------------------------------------%
        beta_i(i,1:end)=atan( ( V(i)+0.5*U_A(i,:) ) ...
            ./( 2*pi*r*n(i)-0.5*U_T(i,:) ));
        d_gamma(i,:)=gammaderive(h*D/2,gamma(i,:));
        
        %--------------------------------------------------------------%
        % Correct the tangential and axial induced velocities with the %
        % induction factor method                                      %
        %--------------------------------------------------------------%
        for j_fixed=2:length(x)-1
            r0=r(j_fixed);
            
            %--------------------------------------------------%
            % Find induction factors over the running variable %
            %--------------------------------------------------%
            for j=1:length(x)
                [i_a(j), i_t(j)] = InductionFactors(r(j),...
                    r0,beta_i(i,j),Z);
            end
            
            %----------------------------------------------------------%
            % Calculate the tangential and axial induced velocities by %
            % integration of Biot-Savarts law with induction factors.  %
            % This requires a singular integration technique           %
            %----------------------------------------------------------%
            U_T(i,j_fixed) = SingularIntegration(r, i_t.*d_gamma(i,:)/(2*pi), r0);
            U_A(i,j_fixed) = SingularIntegration(r, i_a.*d_gamma(i,:)/(2*pi), r0);         
        end
       
        %-----------------------------------------------------------------%
        % Calculate hydrodynamic pitch angle and use it to find effective %
        % angle of attack
        %-----------------------------------------------------------------%
        beta_i(i,1:end)=atan(( V(i)+0.5*U_A(i,:) )./(2*pi.*r*n(i)-0.5*U_T(i,:)) );
        alpha(i,2:end-1) = phi(2:end-1)-beta_i(i,2:end-1);
        
        %-----------------------------------------------------------------%
        % Find lift coeffiecient as the contribution of camber part and   %
        % angle of attack part. Calculate total velocity at foil sections %                                          %
        %-----------------------------------------------------------------%
        C_L(i,:) = C_Lc + 2*pi*alpha(i,:);       
        Vinf(i,:)=sqrt( ( V(i)+0.5*U_A(i,:) ).^2+...
            ( 2*pi*r*n(i)-0.5*U_T(i,:) ).^2);
        
        %------------------------------------------------------%
        % Update the new value of gamma using a damping factor %
        %------------------------------------------------------%
        gamma_new(i,:)=0.5*Vinf(i,:).*chord.*C_L(i,:);
        gamma_input(i,:)=gamma(i,:)+damp*(gamma_new(i,:)-gamma(i,:));
    end
    
    %-----------------------------------------------%
    % Calculate the two-norm error of the iteration %
    %-----------------------------------------------%
    twonorm(iter)= (sum(sum(abs(gamma_input-gamma).^2)))^0.5;
    
    %--------------------------%
    % Updates the distribution %
    %--------------------------%
    gamma=gamma_input;
    
    %---------------------------%
    % Prints percent completion %
    %---------------------------%
    percentdone_print(twonorm(iter)/twonorm(1),tol);
    
    %----------------------------------------------------------------%
    % Checks if convergence is met, if yes then break the iterations %
    %----------------------------------------------------------------%
    if twonorm(iter)/twonorm(1)<tol
        clc
        fprintf('Convergence criterion met after %d iterations. \n',iter);
        break
    end
    
    %------------------------------------------------------------------%
    % Checks if there is any NaN in gamma, if yes then cancels program %
    %------------------------------------------------------------------%
    if isnan(sum(sum(gamma)))==1
        clc
        fprintf('NaN detected! Cancelling program...\n')
        return
    end
end

%-------------------------------------------------------%
% Calculate thrust and torque on every foil section for %
% every advance number                                  %
%-------------------------------------------------------%
for i =1:length(J)
    %---------------------------------------------------------------------%
    % Find thrust and torque contribution because of lift at every foil   %
    % section                                                             %
    %---------------------------------------------------------------------%
    dTi(i,:)=rho*gamma(i,:).*(2*pi*r*n(i)-U_T(i,:)/2);
    dQi(i,:)=rho*gamma(i,:).*(V(i)+U_A(i,:)/2).*r;

    %---------------------------------------------------------------------%
    % Find thrust and torque contribution because of drag at every foil   %
    % section                                                             %
    %---------------------------------------------------------------------%
    dTd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i).*sin(beta_i(i,:));
    dQd(i,:)=0.5*rho*Vinf(i,:).^2.*chord*C_D(i).*r.*cos(beta_i(i,:));

    %--------------------------------------------------------------%
    % Find total thrust and torque component at every foil section %
    %--------------------------------------------------------------%
    dT(i,:)=dTi(i,:)-dTd(i,:);
    dQ(i,:)=dQi(i,:)+dQd(i,:);
end

%--------------------------------------------------%
% Numerical errors at propeller ends are corrected %
%--------------------------------------------------%
dT(:,end)=0;
dQ(:,end)=0;

%--------------------------------------------%
% Integrate with MATLABs trapz() integration %
%--------------------------------------------%
T=Z*trapz(r,dT,2);
Q=Z*trapz(r,dQ,2);

%-------------------------------------------%
% Calculate the thrust and torque coeffient %
%-------------------------------------------%
K_T= T'./(rho*n.^2*D^4);
K_Q= Q'./(rho*n.^2*D^5);

%----------------------%
% Calculate efficiency %
%----------------------%
eta0=J/(2*pi) .* K_T./K_Q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------PLOTTING---------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('rend','painters','pos',[100 100 900 600])
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
