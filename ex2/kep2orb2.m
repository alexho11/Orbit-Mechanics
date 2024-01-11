function [r,nu,GM]=kep2orb(a,e,t,T)
    % kep2orb(a,e,t,T,Omega,omega)
    %
    % Using 6 kepler elements to calculate the polar coordinates of  
    % a satellite and true anomaly.
    % In this case will use time since perigee passage as input,
    % instead of mean anomaly.
    %
    % IN:
    % 6 kepler elements a,e,I,Omega,omega,T
    % a(m) semi-major axis
    % e eccentricity
    % t(s) time since perigee pasasage
    % T(s) time through perigee
    % I(rad) inclination
    % Omega(rad) right ascension of the ascending node
    % omega(rad) argument of perigee
    % T(s) time since perigee passage
    %
    % OUT:
    % r radius(m) 
    % nu true anomaly(rad)               
    % =========================================================
    % author:           Hsin-Feng Ho
    % Martikelnummer:   03770686
    % created at:       26.11.2023
    % last modification:26.11.2023
    % project:          Exercise 1: Keplerian Orbits
    % =========================================================
    
    % gravitational cosntant
    GM=3.986005e14;
    % calculate the mean motion
    n=sqrt(GM./a.^3);
    % calculate the mean anomaly
    M=n.*(t-T);
    % iterate to calculate the eccentric anomaly
    E = M;
    dE = 1;
    for i=1:20
    dE = (M + (e.*sin(E)) - E)./(1 - (e.*cos(E)));
    E = (E + dE);
    end
    E=wrapTo2Pi(E);
    % calculate the true anomaly
    nu=wrapTo2Pi(2*atan2(sqrt(1+e)*sin(E./2), sqrt(1-e)*cos(E./2)));
    % calculate the radius
    r=a.*(1-e.*cos(E));
end