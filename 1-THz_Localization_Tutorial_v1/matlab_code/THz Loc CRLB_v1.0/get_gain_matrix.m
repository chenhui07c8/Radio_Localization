function c = get_gain_matrix(c)
% c.rhoBM = c.GantBM*c.lossBM*c.lambdac/4/pi/c.dBM;
%     c.xiBM = -c.dBM;
    rhoBM_mat = zeros(c.NB, c.NM);
    xiBM_mat = zeros(c.NB, c.NM);
    % complex gain for each channel 
    if c.fc == 60e9
        c.Kabs_mmW = 16;    % 16dB/km
        c.lossBM = 10^(-c.Kabs/1000*c.dBM/10);
    elseif c.fc == 300e9
        c.Kabs_THz = 7.4e-4;% at 3 THz
        c.lossBM = exp(-1/2*c.Kabs_THz*c.dBM);
    end
    if c.gain_type == "Localization"
        % LOS channel
        % Antenna Parameters
        tilde_PhiBM = [c.phiBM_loc c.thetaBM_loc]; % local DOD/DOA from B to R.
        tilde_PhiMB = [c.phiMB_loc c.thetaMB_loc];
        c.GantBM = get_antenna_gain(c.GB, tilde_PhiBM)*get_antenna_gain(c.GM, tilde_PhiMB);
        % channel gain: rho*exp(1j*2*pi/lambda*xi)
        c.rhoBM = c.GantBM*c.lossBM*c.lambdac/4/pi/c.dBM;
        c.xiBM = -c.dBM;
        
        % LOS channel
        if(c.LR > 0)
            % Amplitude and Phase
            if c.fc == 60e9
                c.Kabs_mmW = 16;    % 16dB/km
                c.lossBRM = 10^(-c.Kabs/1000*(c.dRM + c.dBR)/10);
            elseif c.fc == 300e9
                c.Kabs_THz = 7.4e-4;% at 3 THz
                c.lossBRM = exp(-1/2*c.Kabs_THz*(c.dRM + c.dBR));
            end
            % Antenna Parameters
            tilde_PhiBR = [c.phiBR_loc c.thetaBR_loc]; % local DOD/DOA from B to R.
            tilde_PhiMR = [c.phiMR_loc c.thetaMR_loc];
            tilde_PhiRB = [c.phiRB_loc c.thetaRB_loc]; % local DOD/DOA from B to R.
            tilde_PhiRM = [c.phiRM_loc c.thetaRM_loc];
            c.GantBR = get_antenna_gain(c.GB, tilde_PhiBR)*get_antenna_gain(c.GR, tilde_PhiRB);
            c.GantMR = get_antenna_gain(c.GM, tilde_PhiMR)*get_antenna_gain(c.GR, tilde_PhiRM);
            c.GantBRM = c.GantBR*c.GantMR;

            c.rhoBR = c.GantBR*c.lambdac/4/pi/c.dBR;
            c.rhoRM = c.GantMR*c.lambdac/4/pi/c.dRM;
            c.rhoBRM = c.lossBRM*c.rhoBR*c.rhoRM;

            % channel gain: rho*exp(1j*2*pi/lambda*xi)
            c.xiBR = -c.dBR;
            c.xiRM = -c.dRM;
            c.xiBRM = c.xiBR + c.xiRM;
        end
        
        % NLOS channel
        if(c.LC > 0)
            if c.fc == 60e9
                c.Kabs_mmW = 16;    % 16dB/km
                c.lossBCM = 10.^(-c.Kabs_mmW/1000*c.dBCM/10);
            elseif c.fc == 300e9
                c.Kabs_THz = 7.4e-4;% at 3 THz
                c.lossBCM = exp(-1/2*c.Kabs_THz*c.dBCM);
            end
            % geometry
            c.dBC = vecnorm(c.PC-c.PB, 2, 1);
            c.dMC = vecnorm(c.PC-c.PM, 2, 1);
            c.dBCM = c.dBC + c.dMC;
            c.dCM = c.dMC;
            c.tBC = (c.PC - c.PB)./c.dBC;
            [c.phiBC, c.thetaBC] = get_angle_from_dir(c.tBC);
            c.tMC = (c.PC - c.PM)./c.dMC;
            [c.phiMC, c.thetaMC] = get_angle_from_dir(c.tMC);
            c.tCM = -c.tMC;
            [c.phiCM, c.thetaCM] = get_angle_from_dir(c.tCM);

            % BC
            c.tBC_loc = c.RotB'*c.tBC;   % unit direction vector (local) from Tx to Rx
            [c.phiBC_loc, c.thetaBC_loc] = get_angle_from_dir(c.tBC_loc);
            % MC
            c.tMC_loc = c.RotM'*c.tMC;   % unit direction vector (local) from Tx to Rx
            [c.phiMC_loc, c.thetaMC_loc] = get_angle_from_dir(c.tMC_loc);
        
            c.GantBCM = zeros(1, c.LC);
            for lc = 1:c.LC
                tilde_PhiBC = [c.phiBC_loc(lc) c.thetaBC_loc(lc)]; % local DOD/DOA from B to R.
                tilde_PhiMC = [c.phiMC_loc(lc) c.thetaMC_loc(lc)];
                c.GantBCM(lc) = get_antenna_gain(c.GB, tilde_PhiBC)*get_antenna_gain(c.GM, tilde_PhiMC);
            end
            c.rhoBCM = c.lossBCM*c.lambdac/4/pi./(c.dBC + c.dMC).*c.coeClu.*c.GantBCM;
            c.xiBCM = -c.dBCM;
        end
%*******************************************************************************************
    % the rest parts are used for channel comparison.    
    elseif c.AoSA_type == "SSP" || c.AoSA_type == "SPP"
        for b = 1:c.NB
            for m = 1:c.NM
                dbm = norm(c.B(:, b)-c.M(:, m));
                rhoBM_mat(b, m) = c.lambdac/4/pi/dbm;
                xiBM_mat(b, m) = exp(-1j*2*pi/c.lambdac*dbm);
            end
        end
    elseif c.AoSA_type == "PPP"
        for b = 1:c.NB
            for m = 1:c.NM
                dbm = norm(c.PB-c.PM) - (c.M(:,m)-c.PM)'*c.tMB - (c.B(:,b)-c.PB)'*c.tBM;
                rhoBM_mat(b, m) = c.lambdac/4/pi/dbm;
                xiBM_mat(b, m) = exp(-1j*2*pi/c.lambdac*dbm);
            end
        end   
    elseif c.AoSA_type == "SSS"
        for b = 1:c.NB
            for m = 1:c.NM
                gain = get_accurate_channel_element(c, b, m);
                rhoBM_mat(b, m) = abs(gain);
                xiBM_mat(b, m) = gain/abs(gain);
            end
        end
    end
    c.rhoBM_mat = rhoBM_mat;
    c.xiBM_mat = xiBM_mat;
end