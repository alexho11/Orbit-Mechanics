function [rt,rt_dot,az,elev]=efix2topo(r,r_dot)
    % cart2efix(r,r_dot,t)
    %
    % Using cartesian coordinates and time to calculate the cartesian 
    % coordites of a satellite and its velocities in cartesian coordinates
    % in Earth-fixed system
    %
    %
    % IN:
    % cartesian position(s) and velocity(s) of a satellite in inertial frame
    % and its corresponding time.
    %
    % r(m) 3 dimensional vector cartesian coordinates
    % r_dot(m/s) 3 dimensional vector cartesian velocities
    % t(s) corresponding time of the satellite
    % 
    %
    % OUT:
    % re(m) 3-dimensional cartesian coordinates in Earth-fixed frame
    % re_dot(m/s) 3-dimensional cartesian velocities in Earth-fixed frame
    % =========================================================================
    % author:           Hsin-Feng Ho
    % Martikelnummer:   03770686
    % created at:       26.11.2023
    % last modification:26.11.2023
    % project:          Exercise 1: Keplerian Orbits
    % =========================================================================
    
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
    az=atan2(rt(2,:),rt(1,:));
    elev=atan2(rt(3,:),sqrt(rt(1,:).^2+rt(2,:).^2));
end