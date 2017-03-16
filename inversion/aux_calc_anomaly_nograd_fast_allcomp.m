function field = aux_calc_anomaly_nograd_fast_allcomp(PATH_DATA, DATASET, MATRIX, SUSCEPTIBILITIES, OUTPUT_NAME)

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

    BZ_FILENAME = strcat(MODEL_NAME, '_Bz.txt');

    Bz_lines = textread(BZ_FILENAME,'%f', 'commentstyle', 'shell');
    grid_num_pts = length(Bz_lines)/4;

    Bz_vals = reshape(Bz_lines, 4, grid_num_pts)';

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


    %% RESULTING GRIDS - ROWS
    field = sum(MATRIX * SUSCEPTIBILITIES, 2);
    field_bx = field(1:grid_num_pts);
    field_by = field(grid_num_pts+1:grid_num_pts*2);
    field_bz = field(grid_num_pts*2+1:grid_num_pts*3);
   
    res_X_north_cmp = field_bx;
    res_Y_east_cmp = field_by;
    res_Z_vert_cmp = -field_bz;   %%%%NED!!!!!

    sus_string = strrep(mat2str(SUSCEPTIBILITIES'), ' ', '_');
    sus_string = strrep(sus_string, '.', '');
    DATASET = strcat(DATASET);%, '_', sus_string); %%%%%%%%%%%%%%


    %% REUSLTING GRIDS - GRD
    res_X_north_cmp = reshape(res_X_north_cmp,grid_num_lon_pts, grid_num_lat_pts);
    res_Y_east_cmp = reshape(res_Y_east_cmp,grid_num_lon_pts, grid_num_lat_pts);
    res_Z_vert_cmp = reshape(res_Z_vert_cmp,grid_num_lon_pts, grid_num_lat_pts);

    res_X = res_X_north_cmp';
    res_Y = res_Y_east_cmp';
    res_Z = res_Z_vert_cmp';
    

    %% SAVE ALL GRIDS
    lon_coords = linspace(grid_west_border, grid_east_border, grid_num_lon_pts);
    lat_coords = linspace(grid_south_border, grid_north_border, grid_num_lat_pts);

    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bx.xyz'), lon_coords, lat_coords, res_X); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_By.xyz'), lon_coords, lat_coords, res_Y); %%%NEGATIVE
    save_grid_to_xyz_file(strcat(OUTPUT_NAME, '_Bz.xyz'), lon_coords, lat_coords, res_Z); %%%NEGATIVE


end





