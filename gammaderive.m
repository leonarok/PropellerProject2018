function [d_gamma] = gammaderive(h,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----CALCULATES THE DERIVATIVE OF GAMMA BASED ON FINITE DIFFERENCING-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=length(gamma);

%-----------------------------------------------------------------%
% Calculate first and last derivative with backwards and forwards %
% differencing                                                    %
%-----------------------------------------------------------------%
d_gamma(1)=(-gamma(3)+4*gamma(2)-3*...
        gamma(1))/(2*h);
d_gamma(k)=(gamma(k-2)-4*gamma(k-1)+3*...
        gamma(k))/(2*h);

%-------------------------------------------------------%
% Calculate inner derivatives with central differencing %
%-------------------------------------------------------%
for i = 2:k-1
    d_gamma(i)=(-gamma(i-1)+gamma(i+1))/...
        (2*h);
end
             
end
