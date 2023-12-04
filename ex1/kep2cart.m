function [ri,ri_dot]=kep2cart(a,e,t,T,I,Omega,omega)
    % kep2cart(a,e,t,T,I,Omega,omega)
    %
    % Using kepler elements and time to calculate the cartesian coordites of  
    % a satellite and its velocities in cartesian coordinates.
    %
    %
    % IN:
    % Kepler elements a,e,I,Omega,omega
    % a(m) semi-major axis
    % e eccentricity
    % I(rad) inclination
    % Omega(rad) right ascension of the ascending node
    % omega(rad) argument of perigee
    %
    % t(s) time to compute the satellite position
    % T(s) time of perigee
    %
    % OUT:
    % ri(m) 3-dimensional cartesian coordinates in inertial frame
    % ri_dot(m/s) 3-dimensional cartesian velocities in inertial frame
    % =========================================================================
    % author:           Hsin-Feng Ho
    % Martikelnummer:   03770686
    % created at:       26.11.2023
    % last modification:26.11.2023
    % project:          Exercise 1: Keplerian Orbits
    % =========================================================================
     
    % calculate the polar coordinates and true anomaly
    [r,nu,GM]=kep2orb(a,e,t,T);
    
    % calculate the coordinates and velocities in orbit reference frame
    rf=[r.*cos(nu);r.*sin(nu);zeros(1,length(r))];
    rf_dot=[-sqrt(GM./a./(1-e.^2)).*sin(nu);sqrt(GM./a./(1-e.^2)).*(e+cos(nu));zeros(1,length(r))];
    
    % transform to inertial frame
    ri=MatRot(-Omega,3)*MatRot(-I,1)*MatRot(-omega,3)*rf;
    ri_dot=MatRot(-Omega,3)*MatRot(-I,1)*MatRot(-omega,3)*rf_dot;
end