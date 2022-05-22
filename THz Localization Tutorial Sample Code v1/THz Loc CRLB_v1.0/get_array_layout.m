% Generate an ULA or an UPA based on the dimensions
% Sequence of an UPL: left-bottom to left-top
function layout = get_array_layout(N_dim)
    row = N_dim(1); 
    col = N_dim(2);
    [X,Y] = meshgrid((0:(col-1))-(col-1)/2, (0:(row-1))-(row-1)/2);
    layout = [zeros(row*col, 1), X(:), Y(:)]';
end