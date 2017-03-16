function [X, Y, Z] = surface_read_from_file(FILENAME)
    surfData = textread(FILENAME, '%f', 'commentstyle', 'shell');
    nPts = length(surfData)/3;
    Pts = reshape(surfData, 3, nPts)';

    X = Pts(:, 1);
    Y = Pts(:, 2);
    Z = Pts(:, 3);

    return
end

