% convert angle to direction vector or derivative direction vector
function tv = get_dir_from_angle(Phi, type)
    phi = Phi(1,:);
    theta = Phi(2,:);
    
    if(nargin == 1)
        type = 1;
    end

    if type == 1
        tv = [cosd(phi).*cosd(theta); sind(phi).*cosd(theta); sind(theta)];   % steering vector
    elseif type == 2
        tv = [-sind(phi).*cosd(theta); cosd(phi).*cosd(theta); sind(theta)];  % derivative, only phi (azimuth)
    elseif type == 3
        tv = [sind(phi).*sind(theta); -cosd(phi).*sind(theta); cosd(theta)];  % derivative both angles
    end
    
end