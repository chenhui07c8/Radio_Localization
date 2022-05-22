% get the derivative of Phi (phi/theta) from direction vector t
function [D_phi_t, D_theta_t] = get_D_Phi_t(phi, theta)

    D_phi_t = [-sind(phi)*cosd(theta), cosd(phi)*cosd(theta), 0]';
    D_theta_t = [-cosd(phi)*sind(theta), -sind(phi)*sind(theta), cosd(theta)]';

end