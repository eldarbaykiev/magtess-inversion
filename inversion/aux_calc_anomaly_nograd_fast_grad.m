function field = aux_calc_anomaly_nograd_fast_grad(PATH_DATA, DATASET, MATRIX, SUSCEPTIBILITIES, OUTPUT_NAME)

    GMT = 1; %GMT OUTPUT
    %RAD 42588,290808450320743946932468024
    %n = 16 2661, 768
    %n = 15 2839,2193872300213829297954978683

    %% NAME OF DATASET
    PATH = strcat(PATH_DATA, DATASET,'\');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CREATE LIST OF FILES
    FILENAMES = dir(fullfile(strcat(PATH, '*.magtess')));
    N_FILES = length(FILENAMES);


    %% FIRST GRID READING
    MODEL_NAME = strcat(PATH, FILENAMES(1).name);

    BZ_FILENAME = strcat(MODEL_NAME, '_observed_Bxx.xyz');

    Bz_lines = textread(BZ_FILENAME,'%f', 'commentstyle', 'shell');
    grid_num_pts = length(Bz_lines)/3;

    Bz_vals = reshape(Bz_lines, 3, grid_num_pts)';

    grid_west_border = min(Bz_vals(:, 1));
    grid_east_border = max(Bz_vals(:, 1));

    for grid_num_lon_pts = 2 : grid_num_pts
        if(Bz_vals(1, 1) == Bz_vals(grid_num_lon_pts, 1))
            grid_num_lon_pts = grid_num_lon_pts - 1;
            break;
        end
    end

    grid_lon_step = abs(Bz_vals(1, 1) - Bz_vals(2, 1));

    grid_num_lat_pts = grid_num_pts / grid_num_lon_pts;
    grid_south_border = min(Bz_vals(1:grid_num_lat_pts:grid_num_pts, 2));
    grid_north_border = max(Bz_vals(1:grid_num_lat_pts:grid_num_pts, 2));

    grid_lat_step = abs(Bz_vals(1, 2) - Bz_vals(1+grid_num_lon_pts, 2));

%
    %% RESULTING GRIDS - ROWS
    field = sum(MATRIX * SUSCEPTIBILITIES, 2);
    field_bxx = field(1:grid_num_pts);
    field_byx = field(grid_num_pts+1:grid_num_pts*2);
    field_bzx = field(grid_num_pts*2+1:grid_num_pts*3);
    field_bxx = field(grid_num_pts*3+1:grid_num_pts*4);
    field_byx = field(grid_num_pts*4+1:grid_num_pts*5);
    field_bzx = field(grid_num_pts*5+1:grid_num_pts*6);
    
    sus_string = strrep(mat2str(SUSCEPTIBILITIES'), ' ', '_');
    sus_string = strrep(sus_string, '.', '');
    DATASET = strcat(DATASET);%, '_', sus_string); %%%%%%%%%%%%%%


    %% REUSLTING GRIDS - GRD
    res_XX = reshape(field_bxx,grid_num_lon_pts, grid_num_lat_pts);
    res_YX = reshape(field_byx,grid_num_lon_pts, grid_num_lat_pts);
    res_ZX = reshape(field_bzx,grid_num_lon_pts, grid_num_lat_pts);
    res_XY = reshape(field_bxy,grid_num_lon_pts, grid_num_lat_pts);
    res_YY = reshape(field_byy,grid_num_lon_pts, grid_num_lat_pts);
    res_ZY = reshape(field_bzy,grid_num_lon_pts, grid_num_lat_pts);


    %% SAVE ALL GRIDS
    lon_coords = linspace(grid_west_border, grid_east_border, grid_num_lon_pts);
    lat_coords = linspace(grid_south_border, grid_north_border, grid_num_lat_pts);

    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bxx.xyz'), lon_coords, lat_coords, res_XX'); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Byx.xyz'), lon_coords, lat_coords, res_YX'); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bzx.xyz'), lon_coords, lat_coords, res_ZX'); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bxy.xyz'), lon_coords, lat_coords, res_XY'); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Byy.xyz'), lon_coords, lat_coords, res_YY'); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bzy.xyz'), lon_coords, lat_coords, res_ZY'); %%%NEGATIVE


end





