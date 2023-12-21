function [lambda,phi,h]=cart2ell(X,Y,Z,a,e)
%cart2ell(X,Y,Z,a,e)
%
%cell2cart calculates the ellipsoidal coordinates (lambda,phi,h), using the
%ellipsoidal parameters of earth and cartesian coordinates(X,Y,Z,a,e)

% IN:
% X:      x coordinate                                                  [m]
% Y:      y coordinate                                                  [m]
% Z:      z coordinate      
% a:      major axis                                                    [m]
% e:      eccentricity of the ellisoid
%
% OUT:
% lamba:  spherical longtitude                                     [degree]
% phi:    geocentric latitude                                      [degree]
% h:      ellipsoidal height                                            [m]
%
% =========================================================================
% author:           Hsin-Feng Ho
% Martikelnummer:   3378849
% created at:       30.11.2020
% last modification:30.11.2020
% project:          Übung2 Referenzsysteme
% =========================================================================

if nargin ~=5
    error('cart2ell requires 5 arguements');
end
lambda=atan(Y/X)/pi*180;
p=sqrt(X^2+Y^2);

%starting values
phi_0 = atan2(Z,(p*(1-e^2))); % initial geocentric latitude

while true
    N_phi = a/sqrt(1-e^2*sin(phi_0)^2);
    h_i = p/cos(phi_0)-N_phi;
    phi_i = atan2(Z*(N_phi+h_i),p*(N_phi*(1-e^2)+h_i));

    if abs(phi_0-phi_i)<eps
        break;
    end
    phi_0 = phi_i;
end


h=h_i;
phi=phi_i/pi*180