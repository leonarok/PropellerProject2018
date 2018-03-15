function [d_gamma] = gammaderive2(r,gammavec)

P=polyfit(r,gammavec,8);
dPdr=polyder(P);
d_gamma=polyval(dPdr,r);

             
end