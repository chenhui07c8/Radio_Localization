% G is function vector with []
function antennaGain = get_antenna_gain(Gain, Phi)
    antennaGain = 0;
    if length(Phi) == 1     % 2D positioning case
        antennaGain = Gain(1)*(abs(Phi(1))<=Gain(2));
    elseif length(Phi) == 2 % 3D positioning case
        antennaGain = Gain(1)*(abs(Phi(1))<=Gain(2))*(abs(Phi(2))<=Gain(3));
    end
end