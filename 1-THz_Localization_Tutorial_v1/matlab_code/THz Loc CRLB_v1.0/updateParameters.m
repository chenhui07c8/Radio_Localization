function c = updateParameters(c, args)
%     rng(0);
%     c.rng_ind = 6;

    % args: "All", "Geometry", "Channel", "Signal"
    % Geometry: positions, DOA/DOD, directions
    % Channel: AOSA matrix, gain matrix, channel matrix
    % Signal: Tx and Rx symbols

    if(nargin == 1)
        args = "All";
    end

    if c.random_seed == "Random"
        c.rng_ind = randi([0 1000],1);
    end

    % OFDM settings
    c.lambdac = c.lightspeed/c.fc;
    c.dant = c.lambdac/2;    % half wavelength (default antenna spacing)
    c.fdk = -c.BW/2 + c.BW/c.K/2 + (0:(c.K-1))*c.BW/c.K;    % subcarrier frequency
    c.fk = c.fdk + c.fc;
    c.lambdak = c.lightspeed./c.fk;

    % Beam split effect
    if c.beam_split == "True"
        c.beamsplit_coe = c.lambdak/c.lambdac;
    else
        c.beamsplit_coe = ones(size(c.lambdak));
    end

    %********************************************************
    % Update geometry parameters
    %********************************************************
    if args == "All" || args == "Geometry"
    % Check number of channels
        c.NB = prod(c.NB_dim);           % # of BS Elements
        c.NR = prod(c.NR_dim);           % # of RIS Elements
        c.NM = prod(c.NM_dim);           % # of MD Elements

        c.LB = size(c.PB, 2);   % number of BSs (LB LOS paths for each MD)
        c.LM = size(c.PM, 2);   % number of MDs
        c.LC = size(c.PC, 2);   % number of clusters (scatterers, reflectors or diffractors.)
        c.LR = size(c.PR, 2);   % number of RISs

        c.NsaB = prod(c.NsaB_dim);  % number of AE per SA at BS
        c.NsaM = prod(c.NsaM_dim);
        c.NsaR = prod(c.NsaR_dim);

    % Infrastructure & Geometry Parameters
        %************** LOS Geometry ***************
        c.NB = prod(c.NB_dim);      % # of BS Elements (AE for MIMO, SA for AOSA)
        c.NM = prod(c.NM_dim);      % # of MD Elements
        c.NsaB = prod(c.NsaB_dim);  % # of AE per SA at BS
        c.NsaM = prod(c.NsaM_dim);

        % rotation matrix
        c.RotB = eul2rotm(deg2rad(c.OB)', 'ZYX');   % rotation matrix from Euler angles
        c.RotM = eul2rotm(deg2rad(c.OM)', 'ZYX');

        % element local position
        c.B0 = c.DsaB*c.dant*get_array_layout(c.NB_dim);
        c.M0 = c.DsaM*c.dant*get_array_layout(c.NM_dim);

        % element global position
        c.B = c.PB + c.RotB*c.B0;
        c.M = c.PM + c.RotM*c.M0;

        % center distance between BS and MD
        c.dBM = norm(c.PB-c.PM);

        % global DOD/DOA: global = Rotm*local; local = Rotm^-1*global
        c.tBM = (c.PM-c.PB)/c.dBM;
        c.tMB = -c.tBM;
        [c.phiBM, c.thetaBM] = get_angle_from_dir(c.tBM);
        % dir-to-ang:  [atan2d(c.tBM(2,:), c.tBM(1,:)) asind(c.tBM(3,:))]';
        % [cos(c.thetaBM)*cos(c.phiBM); cos(c.thetaBM)*sin(c.phiBM); sin(c.thetaBM)]
        % get_dir_from_angle([c.phiBM; c.thetaBM])

        % local DOD/DOA: global = Rotm*local; local = Rotm^-1*global
        c.tBM_loc = c.RotB'*c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiBM_loc, c.thetaBM_loc] = get_angle_from_dir(c.tBM_loc);
        c.tMB_loc = c.RotM'*-c.tBM;   % unit direction vector (local) from Tx to Rx
        [c.phiMB_loc, c.thetaMB_loc] = get_angle_from_dir(c.tMB_loc);

        %************** RIS Geometry ***************
        if(c.LR > 0)
            c.RotR = eul2rotm(deg2rad(c.OR)', 'ZYX');
            c.R0 = c.DsaR*c.dant*get_array_layout(c.NR_dim);   
            c.R = c.PR + c.RotR*c.R0;
            c.dBR = norm(c.PR-c.PB);
            c.dRB = c.dBR;
            c.dRM = norm(c.PR-c.PM);
            c.dMR = c.dRM;

            % Geometry parameters
            % global DOA/DOD
            c.tBR = (c.PR-c.PB)/c.dBR;
            c.tRB = -c.tBR;
            [c.phiBR, c.thetaBR] = get_angle_from_dir(c.tBR);
            [c.phiRB, c.thetaRB] = get_angle_from_dir(c.tRB);
            c.tRM = (c.PM-c.PR)/c.dRM;
            c.tMR = -c.tRM;
            [c.phiMR, c.thetaMR] = get_angle_from_dir(c.tMR);
            [c.phiRM, c.thetaRM] = get_angle_from_dir(c.tRM);

            % local DOA/DOD
            c.tBR_loc = c.RotB'*c.tBR;   % unit direction vector (local) from Tx to Rx
            [c.phiBR_loc, c.thetaBR_loc] = get_angle_from_dir(c.tBR_loc);
            c.tRB_loc = c.RotR'*c.tRB;   % unit direction vector (local) from Tx to Rx
            [c.phiRB_loc, c.thetaRB_loc] = get_angle_from_dir(c.tRB_loc);
            c.tRM_loc = c.RotR'*c.tRM;   % unit direction vector (local) from Tx to Rx
            [c.phiRM_loc, c.thetaRM_loc] = get_angle_from_dir(c.tRM_loc);
            c.tMR_loc = c.RotM'*c.tMR;   % unit direction vector (local) from Tx to Rx
            [c.phiMR_loc, c.thetaMR_loc] = get_angle_from_dir(c.tMR_loc);
        end

        %************** NLOS Geometry ***************
        if(c.LC > 0)
            % geometry
            c.dBC = vecnorm(c.PC-c.PB, 2, 1);
            c.dMC = vecnorm(c.PC-c.PM, 2, 1);
            c.dBCM = c.dBC + c.dMC;
            c.dCM = c.dMC;
            c.dCB = c.dBC;

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
        end
    end



    %********************************************************
    % Update Channel parameters
    %********************************************************
    if args == "All" || args == "Channel" 
        c = get_Aeq_matrix(c);
        c = get_gain_matrix(c);
        if c.LR > 0
            c.Omega = get_ris_phase(c);   % get ris phase (w/o beamsplit, c.beamsplit)
        end
        c = get_channel(c);
    end


    %********************************************************
    % Update Signal parameters
    %********************************************************

    if args == "All" || args == "Signal"
        % Environment Parameters: False, Sigma, True
        % False: calculate based on Thermal noise: Pn = KTB
        % Sigma: user defined noise level sigma
        % True:  user defined operationBW
        if c.fixed_noise == "False"
            c.operationBW = c.BW/c.K; % True noise should be c.BW, divided by c.K indicates same duration
            c.Pn = c.Kb*c.Temperature*c.operationBW*1000;    % thermal noise linear (in mW)
            c.Pn_dbm = 10*log10(c.Pn);  % thermal noise decibel (in dB)
            c.sigma0 = sqrt(c.Pn);       % Johnson–Nyquist noise: sigma^2 = N_0
            c.NF = 13;   % noise figure 0dB.
            c.sigma = sqrt(10^(c.NF/10))*c.sigma0;      % Output noise level
        elseif c.fixed_noise == "Sigma"
            c.sigma = c.sigma;
        elseif c.fixed_noise == "True"
            c.Pn = c.Kb*c.Temperature*c.operationBW*1000;    % thermal noise linear (in mW)
            c.Pn_dbm = 10*log10(c.Pn);  % thermal noise decibel (in dB)
            c.sigma0 = sqrt(c.Pn);       % Johnson–Nyquist noise: sigma^2 = N_0
            c.NF = 13;   % noise figure 8dB.
            c.sigma = sqrt(10^(c.NF/10))*c.sigma0;      % Output noise level
        end    

        % c.P/c.sigma 
        % c.Omega = c.Omega';
        % size(c.H)

        c = get_Tx_symbols(c);
        c = get_Rx_symbols(c);
    end

end