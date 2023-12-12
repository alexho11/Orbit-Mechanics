function [re,re_dot]=cart2efix(r,r_dot,t)
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
    
    % Earth rotation rate
    Omega_dot=2*pi/86164;
    
    % rotation angle
    theta=Omega_dot.*t+3/12*pi;
    
    re=zeros(size(r,1),size(r,2));
    re_dot=zeros(size(r_dot,1),size(r_dot,2));
    % transformation in Earth fixed system
    for i=1:length(t)
        re(:,i)=MatRot(theta(i),3)*r(:,i);
        re_dot(:,i)=MatRot(theta(i),3)*r_dot(:,i);
    end
    
end