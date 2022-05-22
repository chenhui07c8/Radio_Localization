function J = get_jacobian_matrix_new(p)
    % global: original 8+3 variables
%     D_aM_rhoBM = -p.lambda/4/pi/p.dBM^2*(p.PM-p.PB)/p.dBM;
%     D_aM_thetaBM = [-sin(p.thetaBM)*cos(p.phiBM)/p.dBM; -sin(p.thetaBM)*sin(p.phiBM)/p.dBM; cos(p.thetaBM)/p.dBM];
%     D_aM_phiBM = [-1/p.dBM*sin(p.phiBM)/cos(p.thetaBM); 1/p.dBM*cos(p.phiBM)/cos(p.thetaBM); 0];
%     D_aM_tauBM = 1/p.c*(p.PM-p.PB)/p.dBM;
% 
%     D_aM_rhoBRM = -p.lambda/4/pi/(p.dRM + p.dBR)^2*(p.PM-p.PR)/p.dRM;
%     D_aM_thetaRM = [-sin(p.thetaRM)*cos(p.phiRM)/p.dRM; -sin(p.thetaRM)*sin(p.phiRM)/p.dRM; cos(p.thetaRM)/p.dRM];
%     D_aM_phiRM = [-1/p.dRM*sin(p.phiRM)/cos(p.thetaRM); 1/p.dRM*cos(p.phiRM)/cos(p.thetaRM); 0];
%     D_aM_tauRM = 1/p.c*(p.PM-p.PR)/p.dRM;
% 
%     JPM = [D_aM_rhoBM D_aM_thetaBM D_aM_phiBM D_aM_tauBM D_aM_rhoBRM D_aM_thetaRM D_aM_phiRM D_aM_tauRM];
%     
%     J = zeros(6,11);
%     J(1:3, 1:8) = JPM;
%     J(4:6, 9:11) = eye(3,3);
    

    % global: new 7 variables ((2+2+1)*2 + 1 synchronization + 3 rotation)
    D_aM_rhoBM = -p.lambdac/4/pi/p.dBM^2*(p.PM-p.PB)/p.dBM;
    D_aM_thetaBM = [-sin(p.thetaBM)*cos(p.phiBM)/p.dBM; -sin(p.thetaBM)*sin(p.phiBM)/p.dBM; cos(p.thetaBM)/p.dBM];
    D_aM_phiBM = [-1/p.dBM*sin(p.phiBM)/cos(p.thetaBM); 1/p.dBM*cos(p.phiBM)/cos(p.thetaBM); 0];
    D_aM_dBM = (p.PM-p.PB)/p.dBM;

    D_aM_rhoBRM = -p.lambdac/4/pi/(p.dRM + p.dBR)^2*(p.PM-p.PR)/p.dRM;
    D_aM_thetaRM = [-sin(p.thetaRM)*cos(p.phiRM)/p.dRM; -sin(p.thetaRM)*sin(p.phiRM)/p.dRM; cos(p.thetaRM)/p.dRM];
    D_aM_phiRM = [-1/p.dRM*sin(p.phiRM)/cos(p.thetaRM); 1/p.dRM*cos(p.phiRM)/cos(p.thetaRM); 0];
    D_aM_dRM = (p.PM-p.PR)/p.dRM;

    D_aM_xiBM = -2*pi*p.fc*D_aM_dBM/p.c;
    D_aM_xiBRM = -2*pi*p.fc*D_aM_dRM/p.c;
    

    JPM = [D_aM_rhoBM  D_aM_xiBM  D_aM_thetaBM  D_aM_phiBM  D_aM_dBM...
           D_aM_rhoBRM D_aM_xiBRM D_aM_thetaRM  D_aM_phiRM  D_aM_dRM];

    if p.CRLB_type == "2D4V"   
        % 14: (rho xi theta, phi, tau) 5*2, 3 rotation, 1 clock offset
        % 7: 3 position, 3 rotation, 1 clock offset
        J = zeros(7,14);    
        J(1:3, 1:10) = JPM;
        J(4:6, 11:13) = eye(3,3);
        J(7, 14) = 1;   % syn
        % 4 variables
        J = J([1 2 4 7], [1 2 4 5 6 7 9 10 11 14]);
    elseif p.CRLB_type == "2D8V"
        % 11: 3 position, 3 rotation, 4 gains, 1 clock offset   
        D_aM_rhoBM = zeros(3,1);
        D_aM_rhoBRM = zeros(3,1);
        D_aM_xiBM = zeros(3,1);
        D_aM_xiBRM = zeros(3,1);

        JPM = [D_aM_rhoBM  D_aM_xiBM  D_aM_thetaBM  D_aM_phiBM  D_aM_dBM...
           D_aM_rhoBRM D_aM_xiBRM D_aM_thetaRM  D_aM_phiRM  D_aM_dRM];

        J = zeros(11,14);    
        J(1:3, 1:10) = JPM;
        J(4:6, 11:13) = eye(3,3);
        J(7,1) = 1;
        J(8,2) = 1;
        J(9,6) = 1;
        J(10,7) = 1;
        J(11, 14) = 1;   % syn

        % 8 variables
        J = J([1 2 4 7 8 9 10 11], [1 2 4 5 6 7 9 10 11 14]);
    
    end
    
    
    
end