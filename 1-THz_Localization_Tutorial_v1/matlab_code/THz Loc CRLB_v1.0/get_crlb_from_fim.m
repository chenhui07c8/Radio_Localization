function [CRLB, PEB, OEB, CPEB] = get_crlb_from_fim(FIM, c)
    
    c.J = get_jacobian_matrix(c);
    % size(c.J)
    J = c.J;
    FIM_state = J*FIM*J.';
    if c.syn_type == "Syn"
        CRLB = FIM_state(1:end-1,1:end-1)^(-1);
    else
        CRLB = FIM_state^(-1);
    end


    if c.pos_type + c.ori_type == "2D1D"
        PEB = sqrt(trace(CRLB(1:2,1:2)));
        OEB = sqrt(trace(CRLB(3,3)))/pi*180;
    %     CEB = sqrt(trace(CRLB(end,end)));
    elseif c.pos_type + c.ori_type == "2D0D"
        PEB = sqrt(trace(CRLB(1:2,1:2)));
        OEB = 9999;    
    elseif c.pos_type + c.ori_type == "3D3D"
        PEB = sqrt(trace(CRLB(1:3,1:3)));
        OEB = sqrt(trace(CRLB(4:6,4:6)))/pi*180;
    end
    
    CPEB = zeros(1, c.LC);
    if(c.LC>0)
        for l = 1:c.LC
            if c.pos_type + c.ori_type == "3D3D"
                CPEB(l) = sqrt(trace(CRLB(5*c.LB + 5*c.LR +(1:3)+5*(l-1),5*c.LB + 5*c.LR +(1:3)+5*(l-1))));
            elseif c.pos_type + c.ori_type == "2D1D"
                CPEB(l) = sqrt(trace(CRLB(4*c.LB + 4*c.LR +(1:2)+4*(l-1),4*c.LB + 4*c.LR +(1:2)+4*(l-1))));
            end
%             size(CRLB)
        end
    end
end