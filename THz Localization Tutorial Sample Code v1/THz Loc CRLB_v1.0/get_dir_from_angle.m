% convert angle to direction vector
function tv = get_dir_from_angle(Phi)
    phi = Phi(1);
    theta = Phi(2);
    tv = [cosd(theta)*cosd(phi), cosd(theta)*sind(phi), sind(theta)]';
end