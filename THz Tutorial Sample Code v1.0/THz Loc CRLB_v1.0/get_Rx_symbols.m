% get signal symbols for different channels
% loop K, B, M. (subcarriers, NB, NM)

function c = get_Rx_symbols(c)
% uplink channel
    uBM = zeros(c.NB, c.K);
    uBRM = zeros(c.NB, c.K);
    uBCM = zeros(c.NB, c.K);
    uBCML = zeros(c.NB, c.K, c.LC);

    % LOS
    HL = c.HL;
    for k = 1:c.K
        uBM(:, k) = sqrt(c.P)*HL(:,:,k)*c.Xmk(:,k);
    end
    
    % RIS
    if c.LR > 0
        HR = c.HR;
        for k = 1:c.K
            uBRM(:, k) = sqrt(c.P)*HR(:,:,k)*c.Xmk(:,k);
        end
    end
    
    % NLOS
    if(c.LC > 0)
        HN = c.HN;
        HNL = c.HNL;
        for k = 1:c.K
            uBCM(:,k) = sqrt(c.P)*HN(:,:,k)*c.Xmk(:,k);
            for l = 1:c.LC
                uBCML(:,k,l) = sqrt(c.P)*HNL(:,:,k,l)*c.Xmk(:,k);
            end
        end
    end
    c.uB = uBM + uBRM + uBCM;
    c.uBM = uBM;
    c.uBRM = uBRM;
    c.uBCM = uBCM;
    c.uBCML = uBCML;

    
    % downlink
    uMB = zeros(c.NM, c.K);
    uMRB = zeros(c.NM, c.K);
    uMCB = zeros(c.NM, c.K);
    uMCBL = zeros(c.NM, c.K, c.LC);
    
    for k = 1:c.K
        uMB(:, k) = sqrt(c.P)*transpose(HL(:,:,k))*c.Xbk(:,k);
    end
    
    if c.LR > 0   
        for k = 1:c.K
            uMRB(:, k) = sqrt(c.P)*transpose(HR(:,:,k))*c.Xbk(:,k);
        end
    end

    if(c.LC > 0)
        HN = c.HN;
        HNL = c.HNL;
        for k = 1:c.K
            uMCB(:,k) = sqrt(c.P)*transpose(HN(:,:,k))*c.Xbk(:,k);
            for l = 1:c.LC
                uMCBL(:,k,l) = sqrt(c.P)*transpose(HNL(:,:,k,l))*c.Xbk(:,k);
            end
        end
    end
    c.uM = uMB + uMRB + uMCB;
    c.uMB = uMB;
    c.uMRB = uMRB;
    c.uMCB = uMCB;
    c.uMCBL = uMCBL;
    
end


