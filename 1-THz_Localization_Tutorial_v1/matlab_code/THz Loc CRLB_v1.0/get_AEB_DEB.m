function [CRLB, DEB, AEB] = get_AEB_DEB(FIM, c)
    if c.LB~=1 || c.LR~=0 || c.NM~=1
        error('Double check the number of devices..')
    end
    
    Fm0 = eye(length(FIM));
    if(c.channel_knowledge == "False" && (c.pos_type + c.ori_type == "2D1D"))
        FIM = FIM([1 2 3 5 6 size(FIM,1)], [1 2 3 5 6 size(FIM,1)]);
        if c.syn_type == "Syn"
            CRLB = FIM(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM^(-1);
        end
        AEB = sqrt(trace(CRLB(3,3)))/pi*180;
        DEB = sqrt(trace(CRLB(4,4)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));
    elseif(c.channel_knowledge == "False" && (c.pos_type + c.ori_type == "2D0D"))
        FIM = FIM([1 2 3 5 size(FIM,1)], [1 2 3 5 size(FIM,1)]);
        if c.syn_type == "Syn"
            CRLB = FIM(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM^(-1);
        end
        AEB = sqrt(trace(CRLB(3,3)))/pi*180;
        DEB = sqrt(trace(CRLB(4,4)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));
    elseif((c.channel_knowledge == "Partial") && (c.pos_type + c.ori_type == "2D1D"))
        Fm2 = Fm0([2 3 5 6 size(FIM, 1)], :);
        Fm2(3,1) = -c.lambdac*c.lightspeed/4/pi/c.dBM^2;
%         Fm2(6,5) = -c.lambdac*c.lightspeed/4/pi/(c.dBR + c.dRM)^2;
        FIM2 = Fm2*FIM*Fm2';
        if c.syn_type == "Syn"
            CRLB = FIM2(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM2^(-1);
        end
        AEB = sqrt(trace(CRLB(2,2)))/pi*180;
        DEB = sqrt(trace(CRLB(3,3)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));
    elseif(c.channel_knowledge == "Partial" && (c.pos_type + c.ori_type == "2D0D"))
        Fm2 = Fm0([2 3 5 size(FIM, 1)], :);
        Fm2(3,1) = -c.lambdac*c.lightspeed/4/pi/c.dBM^2;
        FIM2 = Fm2*FIM*Fm2';
        if c.syn_type == "Syn"
            CRLB = FIM2(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM2^(-1);
        end
        AEB = sqrt(trace(CRLB(2,2)))/pi*180;
        DEB = sqrt(trace(CRLB(3,3)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));

    elseif((c.channel_knowledge == "True") && (c.pos_type + c.ori_type == "2D1D"))
        Fm2 = Fm0([3 5 6 size(FIM, 1)], :);
        Fm2(2,1) = -c.lambdac*c.lightspeed/4/pi/c.dBM^2;
        Fm2(2,2) = 2*pi*c.fc;
        FIM2 = Fm2*FIM*Fm2';
        if c.syn_type == "Syn"
            CRLB = FIM2(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM2^(-1);
        end
        AEB = sqrt(trace(CRLB(1,1)))/pi*180;
        DEB = sqrt(trace(CRLB(2,2)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));
     elseif((c.channel_knowledge == "True") && (c.pos_type + c.ori_type == "2D0D"))
        Fm2 = Fm0([3 5 size(FIM, 1)], :);
        Fm2(2,1) = -c.lambdac*c.lightspeed/4/pi/c.dBM^2;
        Fm2(2,2) = 2*pi*c.fc;
        FIM2 = Fm2*FIM*Fm2';
        if c.syn_type == "Syn"
            CRLB = FIM2(1:end-1,1:end-1)^(-1);
        else
            CRLB = FIM2^(-1);
        end
        AEB = sqrt(trace(CRLB(1,1)))/pi*180;
        DEB = sqrt(trace(CRLB(2,2)))*3e8;
%         disp("AEB = " + num2str(AEB) + " | " + "DEB = " + num2str(DEB));
    end
end