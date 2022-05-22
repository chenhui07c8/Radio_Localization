function output = remove_nan(Es, CRBs)

    for i = 1:length(Es)
        if(Es(i) == Inf || Es(i) == -Inf)
            Es(i) = 0.5;
        elseif(isnan(Es(i)))
            Es(i) = CRBs(i);
        end
    end
    output = Es;
end

