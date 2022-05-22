function output = alg_align_index(ind)
    weights = 2.^(1:length(ind));
    ind2 = fliplr(ind);
    while(ind(1)==0)
        ind = circshift(ind, -1);
    end
    while(ind2(1)==0)
        ind2 = circshift(ind2, -1);
    end
    
    if(weights*ind' < weights*ind2')
        output = ind;
    else
        output = ind2;
    end

end