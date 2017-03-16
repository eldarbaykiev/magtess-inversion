function [X, Y, Z] = aux_surface_read_from_file(FILENAME, varargin)

    n_arg = length(varargin);
    if n_arg == 0
        index = 3;     
    else
        index = 4;
    end
 

    
    surfData = textread(FILENAME, '%f', 'commentstyle', 'shell');
    nPts = length(surfData)/index;
    Pts = reshape(surfData, index, nPts)';

    X = Pts(:, 1);
    Y = Pts(:, 2);
    Z = Pts(:, index);

    return
end

