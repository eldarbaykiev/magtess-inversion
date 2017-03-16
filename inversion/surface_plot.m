function surface_plot(Xs, Ys, Zs, varargin)
    
    surf(Xs, Ys, Zs)

    shading flat
    shading interp
    cb = colorbar;
    xlabel(cb, '[m]');
    axis tight;
    ylabel('Y');
    xlabel('X');

    if nargin == 4
        title(varargin{1})
    end
end