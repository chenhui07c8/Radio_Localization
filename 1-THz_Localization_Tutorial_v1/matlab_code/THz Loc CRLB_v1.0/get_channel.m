function c = get_channel(c)
    HL = zeros(c.NB, c.NM, c.K);    % channel matrix
    HR = zeros(c.NB, c.NM, c.K);
    HN = zeros(c.NB, c.NM, c.K);
    HNL = zeros(c.NB, c.NM, c.K, c.LC);

    alphaL = zeros(c.K, 1);         % channel gain to be estimated (localization unknowns)
    XiL = zeros(c.NB, c.NM, c.K);   % extra part of the channel gain (unit complex number)
    
    alphaR = zeros(c.K, 1);
    XiBR = zeros(c.NB, c.NR, c.K);
    XiRM = zeros(c.NR, c.NM, c.K);

    alphaN = zeros(c.K, c.LC);
    XiBC = zeros(c.NB, c.K, c.LC);
    XiMC = zeros(c.NM, c.K, c.LC);
    if c.wave_model == "SWM"
    % LOS channel
        for k = 1:c.K
            for b = 1:c.NB
                for m = 1:c.NM
                    taubm = norm(c.B(:,b) - c.M(:,m))/c.lightspeed;
                    delta_bm = taubm-c.dBM/c.lightspeed;
                    XiL(b, m, k) = exp(-1j*2*pi*(c.fdk(k)*taubm + c.fc*delta_bm));
                end
            end
            alphaL(k) = c.GantBM*c.rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*exp(-1j*2*pi/c.lambdak(k)*c.beta);
            HL(:,:,k) = alphaL(k)*c.AeqBM(:,:,k).*XiL(:,:,k).*c.AeqMB(:,:,k)';
        end
        % RIS channel
        if(c.LR > 0)
            for k = 1:c.K
                for r = 1:c.NR
                    % no AOSA
                    for b = 1:c.NB
                        taubr = norm(c.B(:,b) - c.R(:,r))/c.lightspeed;
                        delta_br = (taubr - c.dBR/c.lightspeed);
                        XiBR(b, r, k) = exp(-1j*2*pi*(c.fdk(k)*taubr + c.fc*delta_br));
                    end
                    for m = 1:c.NM
                        taurm = norm(c.R(:,r) - c.M(:,m))/c.lightspeed;
                        delta_rm = (taurm - c.dRM/c.lightspeed);
                        XiRM(r, m, k) = exp(-1j*2*pi*(c.fdk(k)*taurm + c.fc*delta_rm));
                    end
                end
                alphaR(k) = c.GantBRM*c.rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBRM)*exp(-1j*2*pi/c.lambdak(k)*c.beta);
                for b = 1:c.NB
                    for m = 1:c.NM
                        HR(b, m, k) = alphaR(k)*c.AeqBR(b,:,k).*XiBR(b, :, k)*...
                            diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(XiRM(:, m, k).*c.AeqMR(m, :, k).');
                    end
                end
            end
        end
        % NLOS channel
        if(c.LC > 0)
            for lc = 1:c.LC
                for k = 1:c.K
                    for b = 1:c.NB
                        taubc = norm(c.B(:,b) - c.PC(:,lc))/c.lightspeed;
                        delta_bc = (taubc - c.dBC(lc)/c.lightspeed);
                        XiBC(b, k, lc) = exp(-1j*2*pi*(c.fdk(k)*taubc + c.fc*delta_bc));
                    end
                    for m = 1:c.NM
                        taumc = norm(c.M(:,m) - c.PC(:,lc))/c.lightspeed;
                        delta_mc = (taumc - c.dMC(lc)/c.lightspeed);
                        XiMC(m, k, lc) = exp(-1j*2*pi*(c.fdk(k)*taumc + c.fc*delta_mc));
                    end
                    alphaN(k, lc) = c.GantBCM(lc)*c.rhoBCM(lc)*exp(1j*2*pi/c.lambdac*c.xiBCM(lc))*exp(-1j*2*pi/c.lambdak(k)*c.beta);
                    HNL(:,:,k,lc) = alphaN(k, lc)*(c.AeqBC(:,k,lc).*XiBC(:,k,lc))*(transpose(XiMC(:,k,lc).*c.AeqMC(:,k,lc)));
                end
                HN = HN + HNL(:,:,:,lc);
            end
        end
    elseif c.wave_model == "PWM"
    % LOS channel
        for k = 1:c.K
            for b = 1:c.NB
                for m = 1:c.NM
%                     taubm = norm(c.B(:,b) - c.M(:,m))/c.lightspeed;
                    taubm = (norm(c.PB-c.PM) - (c.M(:,m)-c.PM)'*c.tMB - (c.B(:,b)-c.PB)'*c.tBM)/c.lightspeed;
                    delta_bm = taubm-c.dBM/c.lightspeed;
                    XiL(b, m, k) = exp(-1j*2*pi*(c.fdk(k)*taubm + c.fc*delta_bm));
                end
            end
            alphaL(k) = c.GantBM*c.rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*exp(-1j*2*pi/c.lambdak(k)*c.beta);
            HL(:,:,k) = alphaL(k)*c.AeqBM(:,:,k).*XiL(:,:,k).*c.AeqMB(:,:,k)';
        end
        % RIS channel
        if(c.LR > 0)
            for k = 1:c.K
                for r = 1:c.NR
                    % no AOSA
                    for b = 1:c.NB
%                         taubr = norm(c.B(:,b) - c.R(:,r))/c.lightspeed;
                        taubr = (norm(c.PB-c.PR) - (c.R(:,r)-c.PR)'*c.tRB - (c.B(:,b)-c.PB)'*c.tBR)/c.lightspeed;
                        delta_br = (taubr - c.dBR/c.lightspeed);
                        XiBR(b, r, k) = exp(-1j*2*pi*(c.fdk(k)*taubr + c.fc*delta_br));
                    end
                    for m = 1:c.NM
%                         taurm = norm(c.R(:,r) - c.M(:,m))/c.lightspeed;
                        taurm = (norm(c.PM-c.PR) - (c.R(:,r)-c.PR)'*c.tRM - (c.M(:,m)-c.PM)'*c.tMR)/c.lightspeed;
                        delta_rm = (taurm - c.dRM/c.lightspeed);
                        XiRM(r, m, k) = exp(-1j*2*pi*(c.fdk(k)*taurm + c.fc*delta_rm));
                    end
                end
                alphaR(k) = c.GantBRM*c.rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBRM)*exp(-1j*2*pi/c.lambdak(k)*c.beta);
                for b = 1:c.NB
                    for m = 1:c.NM
                        HR(b, m, k) = alphaR(k)*c.AeqBR(b,:,k).*XiBR(b, :, k)*...
                            diag(c.Omega.*transpose(c.AeqBRM(b, :, m, k)))*(XiRM(:, m, k).*c.AeqMR(m, :, k).');
                    end
                end
            end
        end
        % NLOS channel
        if(c.LC > 0)
            for lc = 1:c.LC
                for k = 1:c.K
                    for b = 1:c.NB
%                         taubc = norm(c.B(:,b) - c.PC(:,lc))/c.lightspeed;
                        taubc = (norm(c.PB-c.PC(:,lc)) - (c.B(:,b)-c.PB)'*c.tBC(:,lc))/c.lightspeed;
                        delta_bc = (taubc - c.dBC(lc)/c.lightspeed);
                        XiBC(b, k, lc) = exp(-1j*2*pi*(c.fdk(k)*taubc + c.fc*delta_bc));
                    end
                    for m = 1:c.NM
%                         taumc = norm(c.M(:,m) - c.PC(:,lc))/c.lightspeed;
                        taumc = (norm(c.PM-c.PC(:,lc)) - (c.M(:,m)-c.PM)'*c.tMC(:,lc))/c.lightspeed;
                        delta_mc = (taumc - c.dMC(lc)/c.lightspeed);
                        XiMC(m, k, lc) = exp(-1j*2*pi*(c.fdk(k)*taumc + c.fc*delta_mc));
                    end
                    alphaN(k, lc) = c.GantBCM(lc)*c.rhoBCM(lc)*exp(1j*2*pi/c.lambdac*c.xiBCM(lc))*exp(-1j*2*pi/c.lambdak(k)*c.beta);
                    HNL(:,:,k,lc) = alphaN(k, lc)*(c.AeqBC(:,k,lc).*XiBC(:,k,lc))*(transpose(XiMC(:,k,lc).*c.AeqMC(:,k,lc)));
                end
                HN = HN + HNL(:,:,:,lc);
            end
        end
    end
    
    c.HL = HL;
    c.alphaL = alphaL;
    c.XiL = XiL;
    
    c.HR = HR;
    c.XiBR = XiBR;
    c.XiRM = XiRM;
    c.alphaR = alphaR;

    c.HN = HN;
    c.HNL = HNL;
    c.XiBC = XiBC;
    c.XiMC = XiMC;
    c.alphaN = alphaN;
    
    c.H = HL + HR + HN;

end