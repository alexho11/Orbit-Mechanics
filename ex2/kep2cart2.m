function [ rc,vc ] = kep2cart2( a,e,r,nu )
GM=398.6005*10^12;
% Computing the position and velocity in the orbital plane-  
rc=[r.*cos(nu);r.*sin(nu);0];
vc=(sqrt(GM/(a*(1-(e^2)))))*[-sin(nu);e+cos(nu);0];
end

