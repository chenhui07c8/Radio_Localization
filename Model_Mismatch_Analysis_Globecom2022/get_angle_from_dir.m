% convert direction to angle
function [phi, theta] = get_angle_from_dir(tv)
    phi = atan2d(tv(2,:), tv(1,:));
    theta = asind(tv(3,:));
end