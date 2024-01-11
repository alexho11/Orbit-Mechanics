function M_R = MatRot(theta, ax)
%MatRot(theta,ax)
%
%MatRot returns a 3*3 rotation matrix of x-/y- or z-axis. The angle and the
%axis are to be given. The unit of the angle is rad. 

%IN:
% theta:     angle of the rotaion                               [rad]
% ax:        1 as x-axis 2 as y-axis 3 as z-axis
%
% OUT:
%M_R         a 3*3 rotation matrix                              
% =========================================================================
% author:           Hsin-Feng Ho
% Martikelnummer:   3378849
% created at:       18.11.2020
% last modification:18.11.2020
% project:          ?bung1 Referenzsysteme
% =========================================================================
switch ax
    case 1 %x-axis rotation matrix
        M_R=[1,0,0;0,cos(theta),sin(theta);0,-sin(theta),cos(theta)];
    case 2 %y-axis rotation matrix
        M_R=[cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
    case 3 %z-axis rotation matrix
        M_R=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
    otherwise
        disp('Error');
end
end
