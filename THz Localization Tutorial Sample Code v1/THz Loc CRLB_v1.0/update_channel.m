function p = update_channel(p, update, arg)
    if exist('update','var')
        for i = 1:length(arg)
            if(arg(i) == "PM")
                p.PM = update{i};
            elseif(arg(i) == "NR")
                p.NR = update{i};
                p.R0 = p.IntR*p.dant*get_array_layout(p.NR, "UPA");
            elseif(arg(i) == "NM")
                p.NM = update{i};
                p.M0 = p.IntM*p.dant*get_array_layout(p.NM, "UPA");
            elseif(arg(i) == "NB")
                p.NB = update{i};
                p.B0 = p.IntB*p.dant*get_array_layout(p.NB, "UPA");
            elseif(arg(i) == "BW")
                p.BW = update{i};
            elseif(arg(i) == "N")
                p.N = update{i};
            end
        end
    end
% 	
%     
%     
    Pn = p.Kb*p.T*p.BW*1000;    % thermal noise linear (in mW)
    Pn_dbm = 10*log10(Pn);  % thermal noise decibel (in dB)
    p.sigma0 = sqrt(Pn);       % Johnson–Nyquist noise: sigma^2 = N_0
    p.NF = 8;   % noise figure 8dB.
    p.sigma = sqrt(10^(p.NF/10))*p.sigma0;       % Johnson–Nyquist noise: sigma^2 = N_0



    p.lambda = p.c/p.fc;
    p.dant = p.lambda/2;    % antenna interval
    p.delta_f = p.BW/p.N;   % subcarrier interval
    p.f0 = p.BW/2;  % central subcarrier
    p.fdn = [1:p.N]*p.BW/p.N;    % subcarrier frequency
    p.fn = p.fdn + p.fc;
    
%     p.fn = p.fn + p.fc;
%     p.f0 = p.f0 + p.fc;

    p.RotB = eul2rotm(p.PhiB, 'ZYX');
    p.RotR = eul2rotm(p.PhiR, 'ZYX');
    p.RotM = eul2rotm(p.PhiM, 'ZYX');

    % element global position
    p.B = p.PB + p.RotB*p.B0;
    p.R = p.PR + p.RotR*p.R0;
    p.M = p.PM + p.RotM*p.M0;

    % center distance
    p.dBM = norm(p.PB-p.PM);
    p.dBR = norm(p.PB-p.PR);
    p.dRM = norm(p.PR-p.PM);

    % global DOD/DOA
    p.tBM = (p.PM-p.PB)/norm(p.PM-p.PB);
    p.phiBM = atan2(p.tBM(2), p.tBM(1));
    p.thetaBM = asin(p.tBM(3));
    % [cos(p.thetaBM)*cos(p.phiBM); cos(p.thetaBM)*sin(p.phiBM); sin(p.thetaBM)]

    p.tRM = (p.PM-p.PR)/norm(p.PM-p.PR);
    p.phiRM = atan2(p.tRM(2), p.tRM(1));
    p.thetaRM = asin(p.tRM(3));
    % [cos(p.thetaRM)*cos(p.phiRM); cos(p.thetaRM)*sin(p.phiRM); sin(p.thetaRM)]

    
    if p.band == "mmWave"
        % d = 1000;
        % pathloss = 10^(-0.016*d/10)
        % 10*log10(pathloss)
        p.Kabs = 16;
        p.lossBM = 10^(-p.Kabs/1000*p.dBM/10);
        p.lossBRM = 10^(-p.Kabs/1000*(p.dRM + p.dBR)/10);
    elseif p.band == "THz"
        p.Kabs = 7.4e-4;
        p.lossBM = exp(-1/2*p.Kabs*p.dBM);
        p.lossBRM = exp(-1/2*p.Kabs*(p.dRM + p.dBR));
    end


    % amplitude
    p.rhoBM = p.KB*p.lossBM*p.lambda/4/pi/p.dBM;
    p.rhoBRM = p.KB*p.KR*p.lossBRM*p.lambda/4/pi/(p.dRM + p.dBR);



    if(size(p.Xmn,1)~=p.N || size(p.Xmn,2)~=p.NM)
        if p.symbol_type == "Ones"
            % fixed symbol
            p.Xmn = ones(p.N, p.NM);   % signal: M antennas, N carriers
            p.Xmn = p.Xmn./abs(p.Xmn);
            p.Xmn = p.Xmn/sqrt(p.NM);
        elseif p.symbol_type == "Random"
            % random x
            rng(0);
            p.Xmn = randn(p.N, p.NM) + 1j*randn(p.N, p.NM);   % signal: M antennas, N carriers
            p.Xmn = p.Xmn./abs(p.Xmn);
            p.Xmn = p.Xmn/sqrt(p.NM);
        end
    end
    p.Omega = get_ris_phase(p);   % get ris phase (w/o beamsplit, p.beamsplit)
%     [uBM, uBRM] = get_received_symbol(p);

end