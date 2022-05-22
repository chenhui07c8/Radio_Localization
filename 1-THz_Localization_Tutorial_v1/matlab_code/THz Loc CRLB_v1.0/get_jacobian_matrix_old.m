function J = get_jacobian_matrix_old(c)

    % global: 5*2 measurements ((2+2+1)*2 + 1 synchronization + 3 rotation)
    D_aM_thetaBM = [-sind(c.thetaBM)*cosd(c.phiBM)/c.dBM; -sind(c.thetaBM)*sind(c.phiBM)/c.dBM; cosd(c.thetaBM)/c.dBM];
    D_aM_phiBM = [-1/c.dBM*sind(c.phiBM)/cosd(c.thetaBM); 1/c.dBM*cosd(c.phiBM)/cosd(c.thetaBM); 0];
    D_aM_tauBM = 1/c.lightspeed*(c.PM-c.PB)/c.dBM;

    D_aM_thetaRM = [-sind(c.thetaRM)*cosd(c.phiRM)/c.dRM; -sind(c.thetaRM)*sind(c.phiRM)/c.dRM; cosd(c.thetaRM)/c.dRM];
    D_aM_phiRM = [-1/c.dRM*sind(c.phiRM)/cosd(c.thetaRM); 1/c.dRM*cosd(c.phiRM)/cosd(c.thetaRM); 0];
    D_aM_tauRM = 1/c.lightspeed*(c.PM-c.PR)/c.dRM;

    D_aM_rhoBM = -c.lambdac/4/pi/c.dBM^2*(c.PM-c.PB)/c.dBM;
    D_aM_rhoBRM = -c.lambdac/4/pi/(c.dRM + c.dBR)^2*(c.PM-c.PR)/c.dRM;

    D_aM_xiBM = -c.lightspeed*D_aM_tauBM;
    D_aM_xiBRM = -c.lightspeed*D_aM_tauRM;
    
    
    %     JPM = [D_aM_rhoBM  D_aM_xiBM  D_aM_thetaBM  D_aM_phiBM  D_aM_tauBM...
    %            D_aM_rhoBRM D_aM_xiBRM D_aM_thetaRM  D_aM_phiRM  D_aM_tauRM];

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
    
    if (c.LC > 0)
%         D_aC_rhoBC = -c.lambdac/4/pi./c.dBC.^2.*(c.PC-c.PB)./c.dBC;
        D_aC_thetaBC = [-sind(c.thetaBC).*cosd(c.phiBC)./c.dBC; -sind(c.thetaBC).*sind(c.phiBC)./c.dBC; cosd(c.thetaBC)./c.dBC];
        D_aC_phiBC = [-1./c.dBC.*sind(c.phiBC)./cosd(c.thetaBC); 1./c.dBC.*cosd(c.phiBC)./cosd(c.thetaBC); zeros(1, c.LC)];
        D_aC_tauBC = 1/c.lightspeed*(c.PC-c.PB)./c.dBC;

%         D_aM_rhoBRM = -c.lambdac/4/pi/(c.dRM + c.dBR)^2*(c.PM-c.PR)/c.dRM;
        D_aM_thetaCM = [-sind(c.thetaCM).*cosd(c.phiCM)./c.dCM; -sind(c.thetaCM).*sind(c.phiCM)./c.dCM; cosd(c.thetaCM)./c.dCM];
        D_aM_phiCM = [-1./c.dCM.*sind(c.phiCM)./cosd(c.thetaCM); 1./c.dCM.*cosd(c.phiCM)./cosd(c.thetaCM); zeros(1, c.LC)];
        D_aM_tauCM = 1/c.lightspeed*(c.PM-c.PC)./c.dCM;
        
        D_aC_thetaCM = [-sind(c.thetaMC).*cosd(c.phiMC)./c.dCM; -sind(c.thetaMC).*sind(c.phiMC)./c.dCM; cosd(c.thetaMC)./c.dCM];
        D_aC_phiCM = [-1./c.dCM.*sind(c.phiMC)./cosd(c.thetaMC); 1./c.dMC.*cosd(c.phiMC)./cosd(c.thetaMC); zeros(1, c.LC)];
        D_aC_tauCM = 1/c.lightspeed*(c.PC-c.PM)./c.dCM;
    end
    
    % unknown channel information, should be most realistic.
    if c.channel_knowledge == "False"
        D_aM_rhoBM = zeros(3,1);
        D_aM_rhoBRM = zeros(3,1);
        D_aM_xiBM = zeros(3,1);
        D_aM_xiBRM = zeros(3,1);
        JPM = [ D_aM_rhoBM  D_aM_xiBM  D_aM_phiBM D_aM_thetaBM D_aM_tauBM...
                D_aM_rhoBRM D_aM_xiBRM D_aM_phiRM D_aM_thetaRM D_aM_tauRM];

        J = zeros(11,14);   % 11 states * 14 measurements    
        J(1:3, 1:10) = JPM;
        J(4:6, 11:13) = eye(3,3);   % rotation
        J(7,1) = 1;     % rhoBM
        J(8,2) = 1;     % xiBM
        J(9,6) = 1;     % rhoBRM
        J(10,7) = 1;    % xiBRM
        J(end, end) = 1;  % syn
        switch c.pos_type + c.ori_type
            case "2D0D"
                if c.LR == 0
                    J = J([1 2 7 8 11], [1 2 3 5 6 7 8 10 11 14]);
                    J(:, 5:9) = 0;
                else
                    J = J([1 2 7 8 9 10 11], [1 2 3 5 6 7 8 10 11 14]);
                    J(:, 9) = 0;
                end
            case "2D1D"
                if c.LR == 0
                % 2 + 1 + 4 + 1 = 8 variables
%                     J = J([1 2 4 7 8 11], [1 2 3 5 size(J,2)-3 size(J,2)]);
                    J = J([1 2 4 7 8 11], [1 2 3 5 6 7 8 10 11 14]);
                    J(:, 5:8) = 0;
                else
                % 2 + 1 + 4 + 1 = 8 variables
                    J = J([1 2 4 7:11], [1 2 3 5 6 7 8 10 11 14]);
%                     J(:, 5:8) = 0;
                end
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 7 8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 7 8 11], :);
                end
            case "3D3D"
                if c.LR == 0
                    J = J([1:8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J;
                end
            otherwise
                error('Position and Orientation formats do not support.');     
        end    
%     
%         if c.pos_type + c.ori_type == "2D1D"
%             % 2 + 1 + 4 + 1 = 8 variables
%             J = J([1 2 4 (7:11)], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
%         elseif c.pos_type + c.ori_type == "3D3D"
%             % 3 + 3 + 4 + 1 =  11 variables
%             J = J;
%         end
        
        if(c.LC > 0)
            J2 = zeros(11 + 5*c.LC, 14 + 8*c.LC);
            for l = 1:c.LC
                JPM_N = [zeros(3, 5) D_aM_phiCM(:,l) D_aM_thetaCM(:,l) D_aM_tauCM(:,l)];
                J2(1:3, (11:18)+8*(l-1)) = JPM_N;
                JPC_N = [zeros(3, 2) D_aC_phiBC(:,l) D_aC_thetaBC(:,l) D_aC_tauBC(:,l) ...
                    D_aC_phiCM(:,l) D_aC_thetaCM(:,l) D_aC_tauCM(:,l)];
                J2((11:13) + 5*(l-1), (11:18)+8*(l-1)) = JPC_N;
                J2(14 + 5*(l-1), 11+8*(l-1)) = 1;   % rhoBCM
                J2(15 + 5*(l-1), 12+8*(l-1)) = 1;   % xiBCM
            end
            
            J2(1:size(JPM,1), 1:size(JPM,2)) = JPM;
            J2(4:6, end-3:end-1) = eye(3,3);
            J2(end, end) = 1;
            % BM + BRM
            J2(7,1) = 1;     % rhoBM
            J2(8,2) = 1;     % xiBM
            J2(9,6) = 1;     % rhoBRM
            J2(10,7) = 1;    % xiBRM
            
            if c.pos_type + c.ori_type == "2D1D"
                % state: 2 + 1 + 4 + [4*c.LC] + 1 variables
                nlos_meas = [11 12 13 15 16 18] + [0:c.LC-1]'*8;
                nlos_meas = reshape(nlos_meas', [1, 6*c.LC]);
                nlos_state = [11 12 14 15] + [0:c.LC-1]'*5;
                nlos_state = reshape(nlos_state', [1, 4*c.LC]);
                J = J2([1 2 4 (7:10) nlos_state size(J2,1)], [1 2 3 5 6 7 8 10 nlos_meas size(J2,2)-3 size(J2,2)]);
            elseif c.pos_type + c.ori_type == "3D3D"
                J = J2;
            end
        end
    elseif c.channel_knowledge == "True"
%         % 14: (rho xi theta, phi, tau) 5*2, 3 rotation, 1 clock offset
%         % 7: 3 position, 3 rotation, 1 clock offset
        JPM = [ D_aM_rhoBM  D_aM_xiBM  D_aM_phiBM D_aM_thetaBM D_aM_tauBM...
                D_aM_rhoBRM D_aM_xiBRM D_aM_phiRM D_aM_thetaRM D_aM_tauRM];

        J = zeros(7,14);    
        J(1:3, 1:10) = JPM;
        J(4:6, 11:13) = eye(3,3);  % rotation
        J(7, 14) = 1;   % syn
%         if c.pos_type + c.ori_type == "2D1D"
%             % 3 + 1 = 4 variables
%             J = J([1 2 4 7], [1 2 3 5 6 7 8 10 11 14]);
%         elseif c.pos_type + c.ori_type == "3D3D"
%             % 6 + 1 = 7 variables
%             J = J;
%         end
        switch c.pos_type + c.ori_type
            case "2D0D"
                if c.LR == 0
                    J = J([1 2 size(J,1)], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                else
                    J = J([1 2 size(J,1)], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                end
            case "2D1D"
                if c.LR == 0
                % 2 + 1 + 4 + 1 = 8 variables
                    J = J([1 2 4 size(J,1)], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                else
                % 2 + 1 + 4 + 1 = 8 variables
                    J = J([1 2 4  size(J,1)], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                end
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 size(J,1)], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 size(J,1)], :);
                end
            case "3D3D"
                if c.LR == 0
                    J = J(:, [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J;
                end
            otherwise
                error('Position and Orientation formats do not support.');     
        end    
    elseif c.channel_knowledge == "Partial"
        % 9: 3 position, 3 rotation, 2 gains, 1 clock offset   
        D_aM_xiBM = zeros(3,1);
        D_aM_xiBRM = zeros(3,1);
        JPM = [ D_aM_rhoBM  D_aM_xiBM  D_aM_phiBM D_aM_thetaBM D_aM_tauBM...
                D_aM_rhoBRM D_aM_xiBRM D_aM_phiRM D_aM_thetaRM D_aM_tauRM];

        J = zeros(11,14);    
        J(1:3, 1:10) = JPM;
        J(4:6, 11:13) = eye(3,3);   % rotation
%         J(7,1) = 1;     % rhoBM
        J(8,2) = 1;     % xiBM
%         J(9,6) = 1;     % rhoBRM
        J(10,7) = 1;    % xiBRM
        J(11, 14) = 1;   % syn
%         if c.pos_type + c.ori_type == "2D1D"
%             % 3 + 4 + 1 = 8 variables
%             J = J([1 2 4 8 10 11], [1 2 3 5 6 7 8 10 11 14]);
%         elseif c.pos_type + c.ori_type == "3D3D"
%             % 3 + 4 + 1 11 variables
%             J = J([1 2 3 4 5 6 8 10 11], :);
%         end
        switch c.pos_type + c.ori_type
            case "2D0D"
                if c.LR == 0
                    J = J([1 2 8 11], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                else
                    J = J([1 2 8 10 11], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                end
            case "2D1D"
                if c.LR == 0
                % 2 + 1 + 4 + 1 = 8 variables
                    J = J([1 2 4 8 11], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                else
                % 2 + 1 + 4 + 1 = 8 variables
                    J = J([1 2 4 8 10 11], [1 2 3 5 6 7 8 10 size(J,2)-3 size(J,2)]);
                end
            case "3D0D"
                if c.LR == 0
                    J = J([1 2 3 7 8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1 2 3 7 8 11], :);
                end
            case "3D3D"
                if c.LR == 0
                    J = J([1:8 11], [1:5 (size(J,2)-3):(size(J,2))]);
                else
                    J = J([1:6 8 10 11], :);
                end
            otherwise
                error('Position and Orientation formats do not support.');     
        end    
    end
    
end