function A = get_beamforming_gain(N, tempN, coe, tilde_Phi, bar_Phi)

    % tempM = pi*deltaM./p.lambdan;
    % with beam squint: coe = 1;
    % without: coe = p.fn/p.fc
    if length(N) == 1 % uniform linear array
        % AN
        
        OmegaN = tempN(1,:).*(sin(tilde_Phi(1))-coe*sin(bar_Phi(1)));    
        AN = 1/sqrt(N(1))*sin(N(1)*OmegaN)./sin(OmegaN);
        AN(OmegaN ==0) = sqrt(N(1));
        AM = 1;
    elseif length(N)==2 % uniform planar array
        % AN
        if tilde_Phi(1) == bar_Phi(1)
            AN = sqrt(N(1));
        else    
            AN = 1/sqrt(N(1))*sin(N(1)*tempN(1,:)*(sin(tilde_Phi(1))*cos(tilde_Phi(2))-sin(tilde_Phi(1))*cos(bar_Phi(1))))...
                /sin(tempN(1,:)*(sin(tilde_Phi(1))*cos(tilde_Phi(2))-sin(tilde_Phi(1))*cos(bar_Phi(1))));
        end
        % AM
        if tilde_Phi(2) == bar_Phi(2)
            AM = sqrt(N(2));
        else    
            AM = 1/sqrt(N(2))*sin(N(2)*tempN(2,:)*(sin(tilde_Phi(2))-sin(bar_Phi(2))));
        end
    end
    
    A = AN.*AM;
end

