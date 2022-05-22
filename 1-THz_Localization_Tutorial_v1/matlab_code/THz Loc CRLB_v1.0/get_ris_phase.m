% core function: optimize RIS phase
% hui.chen@kaust.edu.sa

function Omega = get_ris_phase(c)
    % beamsplit = 0, return N phases.
    % beamsplit = 1, return 1 phase.
    if c.LR > 0
        thetaR = zeros(c.NR, 1);
        if c.RIS_optimizer == "Elzanaty"
            for r = 1:c.NR
                for b = 1:c.NB
                    for m = 1:c.NM
                        taubr = norm(c.B(:,b) - c.R(:,r))/c.lightspeed;
                        taurm = norm(c.R(:,r) - c.M(:,m))/c.lightspeed;
                        taubm = norm(c.B(:,b) - c.M(:,m))/c.lightspeed;
                        thetaR(r) = thetaR(r) + (taubr + taurm - taubm); 
                    end
                end
                %     thetaR(c.R) = thetaR(c.R);
            end
            Omega = exp(1j*thetaR*2*pi.*(c.fc)/c.NM/c.NB);
        elseif c.RIS_optimizer == "He"
            % global: (c.R-c.PR)'*c.tRM || loc: c.R0'*c.tRM_loc
            deltaRM = c.R0'*c.tRM_loc;
            deltaRB = c.R0'*c.tRB_loc;
            Omega = exp(1j*2*pi/c.lambdac*(- deltaRB - deltaRM));
        elseif c.RIS_optimizer == "Random"
            rng(c.rng_ind);
            Omega = exp(-1j*2*pi*rand(c.NR,1));
        elseif c.RIS_optimizer == "Fixed"
            Omega = ones(c.NR,1);
        elseif c.RIS_optimizer == "Quantized"
            for r = 1:c.NR
                for b = 1:c.NB
                    for m = 1:c.NM
                        taubr = norm(c.B(:,b) - c.R(:,r))/c.lightspeed;
                        taurm = norm(c.R(:,r) - c.M(:,m))/c.lightspeed;
                        taubm = norm(c.B(:,b) - c.M(:,m))/c.lightspeed;
                        thetaR(r) = thetaR(r) + (taubr + taurm - taubm); 
                    end
                end
            end
            Omega = exp(1j*thetaR*2*pi.*(c.fc)/c.NM/c.NB);
    %         figure;plot(angle(Omega),'-o');hold on; plot(angle(Omega_discrete),'+');
        elseif c.RIS_optimizer == "Customized"
            Omega = c.Omega;
        end

        if c.RIS_quantization == 1
            partition = 0; % Length 11, to represent 12 intervals
            codebook = (-0.5:1:0.5)*pi; % Length 12, one entry for each interval
            [index, quants] = quantiz(angle(Omega),partition,codebook); % Quantize.
            Omega_discrete = exp(1j*quants);
            Omega = transpose(Omega_discrete);
        elseif c.RIS_quantization == 2
            partition = (-0.5:0.5:0.5)*pi; % Length 11, to represent 12 intervals
            codebook = (-0.75:0.5:0.75)*pi; % Length 12, one entry for each interval
            [index, quants] = quantiz(angle(Omega),partition,codebook); % Quantize.
            Omega_discrete = exp(1j*quants);
            Omega = transpose(Omega_discrete);
        end
    else
        Omega = 0;
    end
end