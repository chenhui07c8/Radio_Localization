% c.multipleTransmission: "Random", "RandomBM", "Fixed"
% "RRR": random beamforming angle at BS/RIS/MD
% "RFR": random beamforming angle at BS/MD, fixed RIS
% % % % "Fixed": fixed beamforming angle at BS/RIS/MD
function cM = get_multiple_transmission_parameters(c)

    rng(c.rng_ind);
    cM.XmkMat = exp(1j*2*pi*rand(c.NM, c.K, c.G));
    cM.XmkMat = cM.XmkMat/norm(cM.XmkMat(:), 'fro');
    % norm(cM.XmkMat(:), 'fro')

    % This part only works for 2D pos.
    if (c.multipleTransmission(1) == 'R')
        cM.BFmatB = zeros(2, c.NB, c.G);
        cM.BFmatB(1,:,:) = (rand(c.NB, c.G)-0.5)*2*90;
    end
    
    if (c.multipleTransmission(3) == 'R')
        cM.BFmatM = zeros(2, c.NM, c.G);
        cM.BFmatM(1,:,:) = (rand(c.NM, c.G)-0.5)*2*90;
    end
    
    if (c.multipleTransmission(2) == 'R') % random
        if(c.LR > 0)
            cM.BFmatRB = zeros(2, c.NR, c.G);
            cM.BFmatRB(1,:,:) = (rand(c.NR, c.G)-0.5)*2*90;
            cM.BFmatRM = zeros(2, c.NR, c.G);
            cM.BFmatRM(1,:,:) = (rand(c.NR, c.G)-0.5)*2*90;
            cM.OmegaMat = exp(-1j*2*pi*rand(c.NR, c.G));
        end
    else % fixed
        if(c.LR > 0)
            cM.BFmatRB = repmat(zeros(2,prod(c.NR_dim)) + [c.phiRB_loc c.thetaRB_loc]', [1, 1, c.G]);
            cM.BFmatRM = repmat(zeros(2,prod(c.NR_dim)) + [c.phiRM_loc c.thetaRM_loc]', [1, 1, c.G]);
            c.RIS_optimizer = "Elzanaty";   % Elzanaty, He, Random, Quantized, Fixed, Customized
            c.Omega = get_ris_phase(c);   % get ris phase (w/o beamsplit, c.beamsplit)
            cM.OmegaMat = repmat(c.Omega, [1, c.G]);
        end
    end
    
end

