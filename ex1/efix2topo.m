function [rt,rt_dot,az,elev]=efix2topo(r,r_dot)
    % efix2topo(r,r_dot) 
    %
    % Using cartesian coordinates and velocities in Earth-fixed 
    % system to calculate the cartesian coordites of a satellite 
    % and its velocities in cartesian coordinates in topocentric 
    % system. In addition to that, the azimuth and elevation of 
    % the satellite are also calculated.
    %
    %
    % IN:
    % cartesian position(s) and velocity(s) of a satellite 
    % in Earth-fixed frame
    % 
    %
    % r(m) 3 dimensional vector cartesian coordinates
    % r_dot(m/s) 3 dimensional vector cartesian velocities
    % 
    %
    % OUT:
    % rt(m) 3-dimensional cartesian coordinates in topocentric frame
    % rt_dot(m/s) 3-dimensional cartesian velocities in topocentric frame
    % az(rad) azimuth of the satellite
    % elev(rad) elevation of the satellite
    % =============================================================
    % author:           Hsin-Feng Ho
    % Martikelnummer:   03770686
    % created at:       02.01.2024
    % last modification:02.01.2024
    % project:          Exercise 1: Keplerian Orbits
    % =============================================================
    
    % translation
    rw=[4075.53022,931.78130,4081.61819]'*1000; % position vector in efix for the station Wettzell
    r_trans=r-rw;
    
    % latitude and longitude of the station Wettzell
    % a=6378137;
    % f_1=298.257223563;
    % f=1/f_1;
    % b=a-a*f;
    % e=sqrt(a^2-b^2)/a;


    % transformation from cartesian to ellipsoidal coordinates
    % [lambda,phi,~]=cart2ell(rw(1),rw(2),rw(3),a,e);
    % lambda=lambda*pi/180;
    % phi=phi*pi/180;
    lambda=12.8781/180*pi;
    phi=49.1449/180*pi;
    % rotation
    Q1=[-1,0,0;0,1,0;0,0,1]; % from right to left-handed system
    rt=Q1*MatRot(pi/2-phi,2)*MatRot(lambda,3)*r_trans;
    rt_dot=Q1*MatRot(pi/2-phi,2)*MatRot(lambda,3)*r_dot;
    az=wrapTo2Pi(atan2(rt(2,:),rt(1,:)));
    elev=atan2(rt(3,:),sqrt(rt(1,:).^2+rt(2,:).^2));
end