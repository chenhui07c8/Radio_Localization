function J = get_jacobian_matrix(c)
    %     JPM = [D_aM_rhoBM  D_aM_xiBM  D_aM_phiBM D_aM_thetaBM D_aM_tauBM...
    %            D_aM_rhoBRM D_aM_xiBRM D_aM_phiRM D_aM_thetaRM D_aM_tauRM];

    % 2D1DPC: 2 Pos + 1 Ori (Syn/Asyn)
    % 2D1DIC: 2 Pos + 1 Ori + 4 (2 complex channel gain) (Syn/Asyn)
    % 3D3DPC: 3 Pos + 3 Ori (Syn/Asyn)
    % 3D3DIC: 3 Pos + 3 Ori + 4 (2 complex channel gain) (Syn/Asyn)

    % measurement vector
    % 3D: 14 variables: 5*(c.LB + c.LR) (rho xi theta phi tau) + 3 ori + 1 syn
    % 2D: 10 variables: 4*(c.LB + c.LR) (rho xi phi tau) + 1 ori + 1 syn
    % 3D with NLOS:     5*c.LB + 5*c.LR + 8*c.LC + 3 + 1
    % 2D with NLOS:     4*c.LB + 6*c.LC + 1 + 1
    
    % state vector
    % 3D: 11 variables: 3 position, 3 rotation, 4 gains, 1 clock offset   
    % 2D: 8 variables:  2 position, 1 rotation, 4 gains, 1 clock offset
    % 3D with NLOS: 10 + c.LC*(3+2) + 1 offset
    % 2D with NLOS: 7 + c.LC*(2+2) + 1 offset
    
    % J_size = length(state)xlength(measure)
    
    % global: 5*2 measurements ((2+2+1)*2 + 1 synchronization + 3 rotation)
    JPM_LOS = zeros(3, 5*c.LB);
    if (c.LB > 0)
        D_aM_thetaBM = [-sind(c.thetaBM)*cosd(c.phiBM)/c.dBM; -sind(c.thetaBM)*sind(c.phiBM)/c.dBM; cosd(c.thetaBM)/c.dBM];
        D_aM_phiBM = [-1/c.dBM*sind(c.phiBM)/cosd(c.thetaBM); 1/c.dBM*cosd(c.phiBM)/cosd(c.thetaBM); 0];
        D_aM_tauBM = 1/c.lightspeed*(c.PM-c.PB)/c.dBM;
        D_aM_rhoBM = -c.lambdac/4/pi/c.dBM^2*(c.PM-c.PB)/c.dBM;
        D_aM_xiBM = -c.lightspeed*D_aM_tauBM;
        for l = 1:c.LB
            JPM_LOS(:, (1:5)+5*(l-1)) = [D_aM_rhoBM(:,l) D_aM_xiBM(:,l)  D_aM_phiBM(:,l) D_aM_thetaBM(:,l) D_aM_tauBM(:,l)];
        end
    end

    JPM_RIS = zeros(3, 5*c.LR);
    if (c.LR > 0)
        D_aM_thetaRM = [-sind(c.thetaRM)*cosd(c.phiRM)/c.dRM; -sind(c.thetaRM)*sind(c.phiRM)/c.dRM; cosd(c.thetaRM)/c.dRM];
        D_aM_phiRM = [-1/c.dRM*sind(c.phiRM)/cosd(c.thetaRM); 1/c.dRM*cosd(c.phiRM)/cosd(c.thetaRM); 0];
        D_aM_tauRM = 1/c.lightspeed*(c.PM-c.PR)/c.dRM;
        D_aM_rhoBRM = -c.lambdac/4/pi/(c.dRM + c.dBR)^2*(c.PM-c.PR)/c.dRM;
        D_aM_xiBRM = -c.lightspeed*D_aM_tauRM;
        for l = 1:c.LR
            JPM_RIS(:, (1:5)+5*(l-1)) = [D_aM_rhoBRM(:,l) D_aM_xiBRM(:,l) D_aM_phiRM(:,l) D_aM_thetaRM(:,l) D_aM_tauRM(:,l)];
        end
    end
    
    JPM_NLOS = zeros(3, 8*c.LC);
    if (c.LC > 0)
        D_aC_thetaBC = [-sind(c.thetaBC).*cosd(c.phiBC)./c.dBC; -sind(c.thetaBC).*sind(c.phiBC)./c.dBC; cosd(c.thetaBC)./c.dBC];
        D_aC_phiBC = [-1./c.dBC.*sind(c.phiBC)./cosd(c.thetaBC); 1./c.dBC.*cosd(c.phiBC)./cosd(c.thetaBC); zeros(1, c.LC)];
        D_aC_tauBC = 1/c.lightspeed*(c.PC-c.PB)./c.dBC;

        D_aM_thetaCM = [-sind(c.thetaCM).*cosd(c.phiCM)./c.dCM; -sind(c.thetaCM).*sind(c.phiCM)./c.dCM; cosd(c.thetaCM)./c.dCM];
        D_aM_phiCM = [-1./c.dCM.*sind(c.phiCM)./cosd(c.thetaCM); 1./c.dCM.*cosd(c.phiCM)./cosd(c.thetaCM); zeros(1, c.LC)];
        D_aM_tauCM = 1/c.lightspeed*(c.PM-c.PC)./c.dCM;
        
        D_aC_thetaCM = [-sind(c.thetaMC).*cosd(c.phiMC)./c.dCM; -sind(c.thetaMC).*sind(c.phiMC)./c.dCM; cosd(c.thetaMC)./c.dCM];
        D_aC_phiCM = [-1./c.dCM.*sind(c.phiMC)./cosd(c.thetaMC); 1./c.dMC.*cosd(c.phiMC)./cosd(c.thetaMC); zeros(1, c.LC)];
        D_aC_tauCM = 1/c.lightspeed*(c.PC-c.PM)./c.dCM;
        
        for l = 1:c.LC
        	JPM_NLOS(:, (1:8) + 8*(l-1)) = [zeros(3, 5) D_aM_phiCM(:,l) D_aM_thetaCM(:,l) D_aM_tauCM(:,l)];
        end
    end
    
    JPM = [JPM_LOS JPM_RIS JPM_NLOS];

    % pos, ori, LOS gain, RIS gain, NLOS pos + NLOS gain
    len_state = 3 + 3 + 2*c.LB + 2*c.LR + 5*c.LC + 1;  
    % 2 gain 2 angle 1 delay (LOS + RIS), 2 gain, phi theta tau (BC), phi theta tau (MC) (NLOS), [3D Near field]
    len_measuremeng = c.LB*5 + c.LR*5 + c.LC*8 + 3 + 1;
    J = zeros(len_state, len_measuremeng);   % 11 states * 14 measurements
    J(4:6, end-3:end-1) = eye(3,3);   % rotation
    J(end, end) = 1;  % syn
    
    % unknown channel information, should be most realistic.
    if c.channel_knowledge == "False"
        % all the gains are treated as unknowns
        zero_gain_ind = [1:5:(5*c.LB + 5*c.LR) 2:5:(5*c.LB + 5*c.LR)...
            (5*c.LB + 5*c.LR + 1):8:(size(J,2)-4) (5*c.LB + 5*c.LR + 2):8:(size(J,2)-4)];
        JPM(:, zero_gain_ind) = 0;
        J(1:3, 1:size(JPM,2)) = JPM;
        for l = 1:c.LB
            J(6 + 2*(l-1) + 1, 5*(l-1) + 1) = 1; % rhoBM
            J(6 + 2*(l-1) + 2, 5*(l-1) + 2) = 1; % xiBM
        end
        for l = 1:c.LR
            J(6 + 2*c.LB + 1 + 2*(l-1), 5*(c.LB + l -1) + 1) = 1; % rhoBRM
            J(6 + 2*c.LB + 2 + 2*(l-1), 5*(c.LB + l -1) + 2) = 1; % xiBRM
        end
        for l = 1:c.LC
            JPC_N = [zeros(3, 2) D_aC_phiBC(:,l) D_aC_thetaBC(:,l) D_aC_tauBC(:,l) ...
                D_aC_phiCM(:,l) D_aC_thetaCM(:,l) D_aC_tauCM(:,l)];
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + (1:3), c.LB*5 + c.LR*5 + 8*(l-1) + (1:8)) = JPC_N;
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 4, c.LB*5 + c.LR*5 + 8*(l-1) + 1) = 1;   % rhoBCM
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 5, c.LB*5 + c.LR*5 + 8*(l-1) + 2) = 1;   % xiBCM
        end
        state_gain_ind = 7:(6+2*(c.LB + c.LR));    % gain for LOS + RIS
        state_nlos_ind = 3 + 3 + 2*c.LB + 2*c.LR + reshape([1 2 4 5]' + 5*(1:c.LC) - 5, [4*c.LC, 1]);
        state_channel_ind = [state_gain_ind state_nlos_ind'];
        switch c.pos_type + c.ori_type
            case "2D0D"
                state_ind = [1 2 state_channel_ind size(J,1)];
                J = J(state_ind, :);
                % set angles (theta) as zero
                zero_angle_ind = [4:5:(5*(c.LB + c.LR)) (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[3 2 1]];
                J(:, zero_ind) = 0;
            case "2D1D"
                % 2D position, 1D angle, LOS_RIS_ind, NLOS_ind, clock
                state_ind = [1 2 4 state_channel_ind size(J,1)];
                J = J(state_ind, :);
                % set angles (theta) as zero
                zero_angle_ind = [4:5:(5*(c.LB + c.LR))...
                    (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4) (5*c.LB + 5*c.LR + 7):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[2 1]];
                J(:, zero_ind) = 0;
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 7 8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 7 8 11], :);
                end
            case "3D3D"
                J = J;
                % set angles (theta) as zero
            otherwise
                error('Position and Orientation formats do not support.');     
        end    

    elseif c.channel_knowledge == "True"
        J(1:3, 1:size(JPM,2)) = JPM;
        for l = 1:c.LC
            JPC_N = [zeros(3, 2) D_aC_phiBC(:,l) D_aC_thetaBC(:,l) D_aC_tauBC(:,l) ...
                D_aC_phiCM(:,l) D_aC_thetaCM(:,l) D_aC_tauCM(:,l)];
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + (1:3), c.LB*5 + c.LR*5 + 8*(l-1) + (1:8)) = JPC_N;
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 4, c.LB*5 + c.LR*5 + 8*(l-1) + 1) = 1;   % rhoBCM
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 5, c.LB*5 + c.LR*5 + 8*(l-1) + 2) = 1;   % xiBCM
        end
        state_gain_ind = [];    % gain for LOS + RIS
        switch c.pos_type + c.ori_type
            case "2D0D"
                state_ind = [1 2 state_gain_ind size(J,1)];
                J = J(state_ind, :);
                zero_angle_ind = [4:5:(5*(c.LB + c.LR)) (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[3 2 1]];
                J(:, zero_ind) = 0;
            case "2D1D"
                % 2D position, 1D angle, LOS_RIS_ind, NLOS_ind, clock
                state_ind = [1 2 4 state_gain_ind size(J,1)];
                J = J(state_ind, :);
                % set angles (theta) as zero
                zero_angle_ind = [4:5:(5*(c.LB + c.LR)) (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[2 1]];
                J(:, zero_ind) = 0;
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 7 8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 7 8 11], :);
                end
            case "3D3D"
                state_ind = [1:6 state_gain_ind size(J,1)];
                if c.LR == 0
                    J = J([1:8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J(state_ind, :);
                end
            otherwise
                error('Position and Orientation formats do not support.');     
        end    
        
    elseif c.channel_knowledge == "Partial"
        zero_gain_ind = 2:5:(5*c.LB + 5*c.LR);  % rho, xi
        JPM(:, zero_gain_ind) = 0;
        J(1:3, 1:size(JPM,2)) = JPM;
        for l = 1:c.LB
            J(6 + 2*(l-1) + 2, 5*(l-1) + 2) = 1; % xiBM
        end
        for l = 1:c.LR
            J(6 + 2*c.LB + 2 + 2*(l-1), 5*(c.LB + l -1) + 2) = 1; % xiBRM
        end
        for l = 1:c.LC
            JPC_N = [zeros(3, 2) D_aC_phiBC(:,l) D_aC_thetaBC(:,l) D_aC_tauBC(:,l) ...
                D_aC_phiCM(:,l) D_aC_thetaCM(:,l) D_aC_tauCM(:,l)];
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + (1:3), c.LB*5 + c.LR*5 + 8*(l-1) + (1:8)) = JPC_N;
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 4, c.LB*5 + c.LR*5 + 8*(l-1) + 1) = 1;   % rhoBCM
            J(3 + 3 + 2*c.LB + 2*c.LR + 5*(l-1) + 5, c.LB*5 + c.LR*5 + 8*(l-1) + 2) = 1;   % xiBCM
        end
        state_gain_ind = 8:2:(6+2*(c.LB + c.LR));    % gain for LOS + RIS
        switch c.pos_type + c.ori_type
            case "2D0D"
                state_ind = [1 2 state_gain_ind size(J,1)];
                J = J(state_ind, :);
                % set angles (theta) as zero
                zero_angle_ind = [4:5:(5*(c.LB + c.LR)) (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[3 2 1]];
                J(:, zero_ind) = 0;
            case "2D1D"
                % 2D position, 1D angle, LOS_RIS_ind, NLOS_ind, clock
                state_ind = [1 2 4 state_gain_ind size(J,1)];
                J = J(state_ind, :);
                % set angles (theta) as zero
                zero_angle_ind = [4:5:(5*(c.LB + c.LR)) (5*c.LB + 5*c.LR + 4):8:(size(J, 2)-4)];
                zero_ind = [zero_angle_ind size(J,2)-[2 1]];
                J(:, zero_ind) = 0;
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 7 8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 7 8 11], :);
                end
            case "3D3D"
                state_ind = [1:6 state_gain_ind size(J,1)];
                if c.LR == 0
                    J = J([1:8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J(state_ind, :);
                end
            otherwise
                error('Position and Orientation formats do not support.');     
        end    

    end
    
end