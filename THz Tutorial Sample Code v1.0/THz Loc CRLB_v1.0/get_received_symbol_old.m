% get signal symbols for different channels
% loop K, B, M. (subcarriers, NB, NM)
function [uB, uBM, uBRM] = get_received_symbol_old(c)
%     uB = zeros()
    uBM = zeros(c.NB, c.K);     % signal in BM channel
    uBRM = zeros(c.NB, c.K);    % signal in uBRM channel
    
    for k = 1:c.K
        for b = 1:c.NB
            ubm = zeros(c.NM, 1);
            ubrm = zeros(c.NM, c.NR);
            for m = 1:c.NM
                taubm = norm(c.B(:,b) - c.M(:,m))/c.lightspeed;
                delta_bm = taubm-c.dBM/c.lightspeed;
                ubm(m) = c.Xmk(m, k)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*exp(-1j*2*pi*(c.fdk(k)*taubm + c.fc*delta_bm));
                % reflect m through NR elements
                for r = 1:c.NR
                    taubr = norm(c.B(:,b) - c.R(:,r))/c.lightspeed;
                    taurm = norm(c.R(:,r) - c.M(:,m))/c.lightspeed;
                    delta_brm = (taubr + taurm - c.dBR/c.lightspeed - c.dRM/c.lightspeed);
                    ubrm(m, r) = c.Xmk(m, k)*exp(-1j*2*pi/c.lambdak(k)*c.beta)*c.Omega(r)*exp(-1j*2*pi*(c.fdk(k)*(taubr + taurm) + c.fc*delta_brm));
                end
            end
            % AOSA
            ubm = ubm.*c.AeqMB(:, b, k).*c.AeqBM(b, :, k)';
            ubrm = ubrm.*c.AeqMR(:, :, k).*c.AeqBR(b, :, k);

            % sum over NM transmitter
            ubBM = c.GantBM*sqrt(c.P)*c.rhoBM*exp(1j*2*pi/c.lambdac*c.xiBM)*sum(ubm);
            ubBRM = c.GantBRM*sqrt(c.P)*c.rhoBRM*exp(1j*2*pi/c.lambdac*c.xiBRM)*sum(ubrm(:));
            
            uBM(b, k) = ubBM;
            uBRM(b, k) = ubBRM;
        end
    end
    uBM = (uBM);
    uBRM = (uBRM);
    uB = uBM + uBRM;
end




% scatterer channel
% uBSM = zeros(p.K, p.NB);
% 
% ubsm = zeros(p.NM, p.NS);
% 
% 
%     for s = 1:p.NS
%         taubs = norm(p.B(:,b) - p.PS(:,s))/p.c;
%         tausm = norm(p.PS(:,s) - p.M(:,m))/p.c;
%         ubsm(m, s) = p.Xmn(n, m)*exp(-1j*2*pi*p.fn(n)*(taubs + tausm));
%     end
% 
% 
% ubBSM = sqrt(p.P)*p.rhoBSM*sum(ubsm(:));
% 
% 
% uBSM(n, b) = ubBSM;
