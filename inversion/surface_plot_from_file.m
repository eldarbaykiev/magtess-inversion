function surface_plot_from_file(FILENAME, varargin)
    [X, Y, Z] = surface_read_from_file(FILENAME);
    
    if nargin == 1
        DISCR = max(abs(min(X)-max(X)), abs(min(Y)-max(Y)));
        for i = 1 : length(Z)-1
            for j = i + 1 : length(Z)
                x1 = X(i);
                y1 = Y(i);
                x2 = X(j);
                y2 = Y(j);
                dist = sqrt((x1-x2)^2 + (y1-y2)^2);
                if dist < DISCR
                    DISCR = dist;
                end
            end
        end
        DISCR = DISCR/2;
    else
        DISCR = varargin{1};
    end

    F = scatteredInterpolant(X,Y,Z);
    [Xs, Ys] = meshgrid(min(X):DISCR:max(X), min(Y):DISCR:max(Y));
    Zs = F(Xs, Ys);

    surface_plot(Xs, Ys, Zs, FILENAME)
end