function c = get_Tx_symbols(c)

    if c.symbol_type == "Ones"
        % fixed symbol
        c.Xmk = ones(c.NM, c.K);   
        c.Xmk = c.Xmk./abs(c.Xmk);
        c.Xmk = c.Xmk/sqrt(c.NM)/sqrt(c.K);
        c.Xbk = ones(c.NB, c.K);   
        c.Xbk = c.Xbk./abs(c.Xbk);
        c.Xbk = c.Xbk/sqrt(c.NB)/sqrt(c.K);
    % norm(c.Xmk)
    elseif c.symbol_type == "Random"
        rng(c.rng_ind);
        c.Xmk = exp(1j*2*pi*rand(c.NM, c.K));
        c.Xmk = c.Xmk/sqrt(c.NM)/sqrt(c.K);
        c.Xbk = exp(1j*2*pi*rand(c.NB, c.K));   % signal: M antennas, N carriers
        c.Xbk = c.Xbk/sqrt(c.NB)/sqrt(c.K);
%       figure;plot(c.Xmk,'o')
%       norm(c.Xmk, 'Fro')
%       sum(sum(abs(c.Xmk(:))))
    elseif c.symbol_type == "RandomSymbol"
        % random symbols but constant beamforming
        w = ones(c.NM,1);
        symbs = exp(1j*2*pi*rand(1,c.K));
        c.Xmk = w*symbs;
        c.Xmk = c.Xmk./abs(c.Xmk);
        c.Xmk = c.Xmk/sqrt(c.NM)/sqrt(c.K);
        w = ones(c.NB,1);
        symbs = exp(1j*2*pi*rand(1,c.K));
        c.Xbk = w*symbs;
        c.Xbk = c.Xbk./abs(c.Xbk);
        c.Xbk = c.Xbk/sqrt(c.NB)/sqrt(c.K);        
    elseif c.symbol_type == "Customized" 
        c.Xmk = c.Xmk;
        c.Xbk = c.Xbk;
    end
    
    
end