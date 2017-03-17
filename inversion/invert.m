%% MAGNETIC TESSEROID INVERSION ALGORITHM
%Eldar Baykiev, 2017

%This script uses IGRF calculator from https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model

%clear everything before start
clc;
clear all;
close all;

%Constants
EARTH_RADIUS_IGRF_KM = 6371.2;
EARTH_RADIUS_TESS_M = 6378137.0;

%STEPS TO RUN
STEP1 = 1; %1 - on; 0 - off  %creation of tesseroid model
STEP2 = 1; %1 - on; 0 - off  %calculation of each tesseroid effect
STEP3 = 1; %1 - on; 0 - off  %design matrix creation
STEP4 = 1; %1 - on; 0 - off  %vector field comnponents inversion
STEP5 = 1; %1 - on; 0 - off  %gradient components inversion

%block filename template
BLOCK_FILENAME_TEMPLATE = 'D%.1f_W%.1fE%.1fS%.1fN%.1f_%s_%s_layer_%d_%d.magtess_block';

%Current date
CURRENT_DATE = datestr(clock, 'ddmmyy_hhMMSS');

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IRGF11 date   
%Specify date that would be used in IGRF11 to calculate the inducing
%magnetic field
IGRF_DATE = '01-Jan-2014';

%Path to surface files
%Specify path to the surface files, that would be used to define the 
%geometry of tesseroid model. Path must have \ in the end.
PATH_SURF = 'C:\data_folder\';

%Top surface filename (in the folder PATH). Format: lon [deg] lat [deg]
%depth [m]. Use a simple name without -, + symbols. This file should be in
%PATH_SURF folder
TOP_FILENAME = 'Surface_Top.xyz';

%Bottom surface filename (in the folder PATH). Format: lon [deg] lat [deg]
%depth [m]. Use a simple name without -, + symbols. This file should be in
%PATH_SURF folder
BOT_FILENAME = 'Surface_Bot.xyz';

%Path to the output
%Specify folder where all intermediate files would be saved
PATH_OUTPUT = 'C:\INVERSION_RESULTS\';

%Tesseroid width (in [deg])
%Set what width tesseroids in the model should have
TESSEROID_WIDTH = 2;

%Tesseroid layer height (in [m])
%Set tesseroid layer's height. Model would be divided into these layers
LAYER_HEIGHT = 20000;

%main edges (in geographical coord)
WEST_EDGE = -180;
EAST_EDGE = 180;
SOUTH_EDGE = -90;
NORTH_EDGE = -55;

%edge extension (in [deg]). All should be positive
EDGE_EXTENSION_W = 0;
EDGE_EXTENSION_E = 0;
EDGE_EXTENSION_S = 0;
EDGE_EXTENSION_N = 5;

%altitude of the resulting grid(in [m])
ALTITUDE = 400000; 

%resulting grids' spacing (in [deg])
SPACING = 2;

%observed data
%Bz in NED coordinate system. Use forward_calc_glob_grid_from_gh_shc
OBSERVED_DATA_FILENAME_BX = strcat('C:\data_folder\observed_Bx.txt');
OBSERVED_DATA_FILENAME_BY = strcat('C:\data_folder\observed_By.txt');
OBSERVED_DATA_FILENAME_BZ = strcat('C:\data_folder\observed_Bz.txt'); %Bz in NED coordinate system. 

%Vector comp. inverison result folder. Result folder would in PATH_OUTPUT
%folder
RESULT_FOLDER_VECT = strcat('VECT_INV\');

%Gradient inverison result folder. Result folder would in PATH_OUTPUT
%folder
RESULT_FOLDER_GRAD = strcat('GRAD_INV\');

%a-priori susceptibility of layers
%note that layers appear in order, that is shown in matlab with command
%"dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));" (see around line 68)
STANDART_SUSCEPTIBILITY = [0.02, 0.02, 0.02, 0.02];

%here you must specify a-priori susceptibility grids for corresponding
%layers. '%d' should be in place where the layer number is in the filename
APRIORI_SUSCEPT_FILENAME_TEMPLATE = strcat('INVERSION_CODE\\data\\suscept_layer%d.xyz');

%layers' susceptibility variances
sigma_x = [0.01; 0.01; 0.01; 0.01]; %layers

%observed data variances
sigma_vect_d = [0.5; 0.5; 0.5]; %Bx By Bz nT

%observed data variances
sigma_grad_d = [5; 5; 5; 5; 5; 5;]; %Bxx Byx Bzx Bxy Byy Bzy pT/km


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist(strcat(PATH_SURF)))
    fprintf('Folder with surfaces %s does not exist. \n', PATH_SURF);
    return;
end

if (PATH_SURF(end) ~= '\')
    PATH_SURF = strcat(PATH_SURF, '\');
    fprintf('Missing \\ added to the end of PATH_SURF. \n');
end

if (~exist(strcat(PATH_SURF, TOP_FILENAME)))
    fprintf('Top surface file %s is not in the folder %s. \n', TOP_FILENAME, PATH_SURF);
    return;
end
    
if (~exist(strcat(PATH_SURF, BOT_FILENAME)))
    fprintf('Top surface file %s is not in the folder %s. \n', BOT_FILENAME, PATH_SURF);
    return;
end

if (~exist(PATH_OUTPUT))
    mkdir(PATH_OUTPUT)
end

if (PATH_OUTPUT(end) ~= '\')
    PATH_OUTPUT = strcat(PATH_SURF, '\');
    fprintf('Missing \\ added to the end of PATH_OUTPUT. \n');
end


%% %%%%%%%%%%%%%%%%%   STEP 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This step uses two surfaces and IGRF11 to create magnetized multilayered tesseroid model
if STEP1
    close all;
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    
    fprintf('STEP 1\n');
    %% LOADING SURFACES
    
    figure(1);
    surface_plot_from_file(strcat(PATH_SURF, TOP_FILENAME), 1);
    hold on
    surface_plot_from_file(strcat(PATH_SURF, BOT_FILENAME), 1);
    shg

    [X, Y, Ztop] = surface_read_from_file(strcat(PATH_SURF, TOP_FILENAME));
    F_Ztop = scatteredInterpolant(X, Y, Ztop);

    [X, Y, Zbot] = surface_read_from_file(strcat(PATH_SURF, BOT_FILENAME));
    F_Zbot = scatteredInterpolant(X, Y, Zbot);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    IGRF_DATENUM = datenum(IGRF_DATE);

    B_W = WEST_EDGE -EDGE_EXTENSION_W; 
    B_E = EAST_EDGE +EDGE_EXTENSION_E; 
    B_S = SOUTH_EDGE-EDGE_EXTENSION_S; 
    B_N = NORTH_EDGE+EDGE_EXTENSION_N; 

    MAX_BOTTOM_DEPTH = floor(min(Zbot)/10000.0)*10000.0;
    MAX_TOP_DEPTH = ceil(max(Ztop)/10000.0)*10000.0;


    for range_of_depth = MAX_BOTTOM_DEPTH : LAYER_HEIGHT : MAX_TOP_DEPTH-LAYER_HEIGHT/2.0;

        for i = 1 : length(range_of_depth)

            HOB = range_of_depth(i);
            HOT = range_of_depth(i) + LAYER_HEIGHT; 

            TESS_FILENAME = sprintf(BLOCK_FILENAME_TEMPLATE, TESSEROID_WIDTH, B_W, B_E, B_S, B_N, TOP_FILENAME, BOT_FILENAME, HOB, HOT);
            TESS_FILE = fopen(strcat(PATH_OUTPUT, TESS_FILENAME), 'w');

            for blon = B_W : TESSEROID_WIDTH : B_E-TESSEROID_WIDTH/2.0
                for blat = B_S : TESSEROID_WIDTH : B_N-TESSEROID_WIDTH/2.0
                    HOB = range_of_depth(i);
                    HOT = range_of_depth(i) + LAYER_HEIGHT; 


                    W = blon;
                    E = blon + TESSEROID_WIDTH;
                    S = blat;
                    N = blat + TESSEROID_WIDTH;

                    lon_c = 0.5*(W + E);
                    lat_c = 0.5*(S + N);

                    curr_Zbot = F_Zbot(lon_c, lat_c);
                    curr_Ztop = F_Ztop(lon_c, lat_c);

                    if(curr_Zbot > HOB)
                        HOB = curr_Zbot;
                    end

                    if(curr_Ztop < HOT)
                        HOT = curr_Ztop;
                    end

                    if HOT <= HOB
                        continue
                    else
                        alt_c = 0.5*(HOT + HOB);

                        altitude_igrf = (alt_c + EARTH_RADIUS_TESS_M - EARTH_RADIUS_IGRF_KM*1000.0)/1000.0;           
                        [Bx, By, Bz] = igrf(IGRF_DATENUM, lat_c, lon_c, altitude_igrf);
                        Bz = - Bz;

                        fprintf(TESS_FILE, '%f %f %f %f %f %f %d %d %.10f %.10f %.10f\n', W, E, S, N, HOT, HOB, i, 1, Bx, By, Bz);  
                    end

                end          
            end

            fclose(TESS_FILE);
        end

    end

    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 1 finished\n');
end


%% %%%%%%%%%%%%%%%%%   STEP 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script calculates the magnetic effect of each tesseroid in blocks
% note: before running this script, put tessbx, tessby and
% tessbz into the output folder PATH_OUTPUT

if STEP2
    close all;
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    
    fprintf('STEP 2\n');
    
    if( ~(exist(strcat(PATH_OUTPUT, 'tessbx.exe')) && exist(strcat(PATH_OUTPUT, 'tessby.exe')) && exist(strcat(PATH_OUTPUT, 'tessbz.exe')) && exist(strcat(PATH_OUTPUT, 'tessgrd.exe')) ))
    	fprintf('Missing tessbx.exe, tessby.exe, tessbz.exe or tessgrd.exe. Please copy these files to the output directory %s \n', PATH_OUTPUT);
        return
    end

    %% LOADING BLOCKS
    
    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);

    fprintf('Files to process:\n');
    for n = 1 : N_FILES
        fprintf('  %d. %s\n', n, FILENAMES(n).name);
    end

    %% CALCULATING
    fprintf('Calculating...\n')


    for i = 1 : N_FILES
        WEST = WEST_EDGE;              
        EAST = EAST_EDGE;              
        SOUTH = SOUTH_EDGE;            
        NORTH = NORTH_EDGE; 
        
        current = FILENAMES(i).name;
        BLOCK = sprintf('%s\\', current(1:end-14));
        mkdir(strcat(PATH_OUTPUT, BLOCK));

        N_LON = length(WEST : SPACING : EAST);             
        N_LAT = length(SOUTH : SPACING : NORTH);  

        fprintf('Creating grid (w:%f e:%f s:%f n:%f nlon:%f nlat:%f alt:%f)... ', WEST, EAST, SOUTH, NORTH, N_LON, N_LAT, ALTITUDE);
        system(strcat(PATH_OUTPUT, 'tessgrd.exe -r', num2str(WEST), '/', num2str(EAST), '/', num2str(SOUTH), '/', num2str(NORTH), ' -b', num2str(N_LON), '/', num2str(N_LAT), ' -z', num2str(ALTITUDE), ' > ', PATH_OUTPUT, BLOCK, 'grid_earth.txt'));
        fprintf('finished\n');

        fprintf('  %d. %s \n', i, current);


        Nrows = numel(textread(strcat(PATH_OUTPUT, current),'%1c%*[^\n]'));
        h = waitbar(0,sprintf('Calculating file no. %d...', i));

        fid = fopen(strcat(PATH_OUTPUT, current), 'r');

        j = 1;
        tline = fgetl(fid);
        while ischar(tline)
            perc = (j/Nrows);
            waitbar(perc,h,sprintf('Calculating file no. %d...\n%.2f%%', i,perc*100));


            fid_temp = fopen(strcat(PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess'), 'w');
            fprintf(fid_temp, '%s', tline);
            fclose(fid_temp);

            command = '';
            command_bx = strcat(PATH_OUTPUT, 'tessbx.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_Bx.txt');
            command_by = strcat(PATH_OUTPUT, 'tessby.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_By.txt');
            command_bz = strcat(PATH_OUTPUT, 'tessbz.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_Bz.txt');   

            tline = fgetl(fid);
            j = j+1;
            if ischar(tline)
                fid_temp = fopen(strcat(PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess'), 'w');
                fprintf(fid_temp, '%s', tline);
                fclose(fid_temp);

                command_bx2 = strcat(PATH_OUTPUT, 'tessbx.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_Bx.txt');
                command_by2 = strcat(PATH_OUTPUT, 'tessby.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_By.txt');
                command_bz2 = strcat(PATH_OUTPUT, 'tessbz.exe "', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess', '" < ', PATH_OUTPUT, BLOCK, 'grid_earth.txt > ', PATH_OUTPUT, BLOCK, 'tess_n', num2str(j), '.magtess_Bz.txt');   

                command = strcat(command_bx, ' | ', command_by, ' | ',command_bz, ' | ',command_bx2, ' | ', command_by2, ' | ',command_bz2);
                tline = fgetl(fid);
                j = j+1;
            else
                command = strcat(command_bx, ' | ', command_by, ' | ',command_bz);
            end


            system(command);
        end

        close(h);
        fclose(fid);



    end
    
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 2 finished\n');
end

%% %%%%%%%%%%%%%%%%%   STEP 3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if STEP3
    close all;
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 3\n');
    
    %% This script calculates all design matrices
    
    
    %% LOADING BLOCKS
    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);

    fprintf('Blocks to process:\n');
    for n = 1 : N_FILES
        fprintf('  %d. %s\n', n, FILENAMES(n).name);
    end

    %% CREATE GRADIENT GRIDS
    
    WEST = WEST_EDGE;              
    EAST = EAST_EDGE;              
    SOUTH = SOUTH_EDGE;            
    NORTH = NORTH_EDGE; 
        
    N_LON = length(WEST : SPACING : EAST);             
    N_LAT = length(SOUTH : SPACING : NORTH);  

    TEST_SUM_Bxx = zeros(N_LAT-2, N_LON-2);
    TEST_SUM_Byx = zeros(N_LAT-2, N_LON-2);
    TEST_SUM_Bzx = zeros(N_LAT-2, N_LON-2);
    TEST_SUM_Bxy = zeros(N_LAT-2, N_LON-2);
    TEST_SUM_Byy = zeros(N_LAT-2, N_LON-2);
    TEST_SUM_Bzy = zeros(N_LAT-2, N_LON-2);



    %% CALCULATING
    for hhh = 1 : N_FILES
        %recognize corresponding folder


        current = FILENAMES(hhh).name;
        DATASET = sprintf('%s', current(1:end-14));
        sprintf('Current folder: %s\n',DATASET);

        PATH_OUTPUT_BLOCK = strcat(PATH_OUTPUT, DATASET,'\');


        BODIES = dir(fullfile(strcat(PATH_OUTPUT_BLOCK, '*.magtess')));
        N_BODIES = length(BODIES);
        if(N_BODIES == 0)
            printf('No calculated effects of individual tesseroids\n')
            return
        end

        dummy = textread(strcat(PATH_OUTPUT_BLOCK, BODIES(1).name, '_Bz.txt'),'%f', 'commentstyle', 'shell');
        grid_num_pts = length(dummy)/4;
        % MAGNETIC FIELD COMP

        A = zeros(grid_num_pts*3, N_BODIES);
        BAR_TEXT = 'Reading design matrix...';

        h = waitbar(0,'BAR_TEXT');
        for i = 1 : N_BODIES
            %PROGRESS BAR
            if rem(i, 100) == 0
                perc = (i/N_BODIES);
                waitbar(perc,h,sprintf(strcat(BAR_TEXT, '\n%.1f%%'),perc*100));
            end


            [BX_VALS, BY_VALS, BZ_VALS] = load_vals(strcat(PATH_OUTPUT_BLOCK, BODIES(i).name));
            A(:, i) = [BX_VALS(:); BY_VALS(:); BZ_VALS(:)];

        end
        waitbar(1,h,sprintf(strcat(BAR_TEXT, '\n100%%')));
        close(h);

        fileID = fopen(strcat(PATH_OUTPUT_BLOCK, 'design_matrix.bin'),'w');
        fwrite(fileID,A,'float32');
        fclose(fileID); 

        % GRADIENT

        dummy_VALS = reshape(dummy, 4, grid_num_pts)';

        grid_west_border = min(dummy_VALS(:, 1));
        grid_east_border = max(dummy_VALS(:, 1));

        for grid_num_lon_pts = 2 : grid_num_pts
            if(dummy_VALS(1, 1) == dummy_VALS(grid_num_lon_pts, 1))
                grid_num_lon_pts = grid_num_lon_pts - 1;
                break;
            end
        end

        grid_lon_step = abs(dummy_VALS(1, 1) - dummy_VALS(2, 1));

        grid_num_lat_pts = grid_num_pts / grid_num_lon_pts;
        grid_south_border = min(dummy_VALS(1:grid_num_lat_pts:grid_num_pts, 2));
        grid_north_border = max(dummy_VALS(1:grid_num_lat_pts:grid_num_pts, 2));

        grid_lat_step = abs(dummy_VALS(1, 2) - dummy_VALS(1+grid_num_lon_pts, 2));

        grid_height = dummy_VALS(1, 3);

        A = zeros((grid_num_lat_pts-2)*(grid_num_lon_pts-2)*6, N_BODIES);
        BAR_TEXT = 'Reading design matrix...';

        h = waitbar(0,'BAR_TEXT');
        for ni = 1 : N_BODIES
            %PROGRESS BAR
            if rem(ni, 100) == 0
                perc = (ni/N_BODIES);
                waitbar(perc,h,sprintf(strcat(BAR_TEXT, '\n%.1f%%'),perc*100));
            end

            OBSERVED_NAME = strcat(PATH_OUTPUT_BLOCK, BODIES(ni).name);
            OBSERVED_BX_FILENAME = strcat(OBSERVED_NAME, '_Bx.txt');
            OBSERVED_BY_FILENAME = strcat(OBSERVED_NAME, '_By.txt');
            OBSERVED_BZ_FILENAME = strcat(OBSERVED_NAME, '_Bz.txt');

            OBSERVED_BX = textread(OBSERVED_BX_FILENAME,'%f', 'commentstyle', 'shell');
            OBSERVED_BY = textread(OBSERVED_BY_FILENAME,'%f', 'commentstyle', 'shell');
            OBSERVED_BZ = textread(OBSERVED_BZ_FILENAME,'%f', 'commentstyle', 'shell');
            grid_num_pts = length(OBSERVED_BX)/4;

            BX_VALS = reshape(OBSERVED_BX, 4, grid_num_pts)';
            BY_VALS = reshape(OBSERVED_BY, 4, grid_num_pts)';
            BZ_VALS = reshape(OBSERVED_BZ, 4, grid_num_pts)';

            grid_field_bx = BX_VALS(:, 4);
            grid_field_by = BY_VALS(:, 4);
            grid_field_bz = BZ_VALS(:, 4); 

            field_bx = zeros(grid_num_pts, 1);
            field_by = zeros(grid_num_pts, 1);
            field_bz = zeros(grid_num_pts, 1);

            field_bx = field_bx + grid_field_bx;
            field_by = field_by + grid_field_by;
            field_bz = field_bz + grid_field_bz; 

            res_X_north_cmp = field_bx;
            res_Y_east_cmp = field_by;
            res_Z_vert_cmp = -field_bz;

            res_X_north_cmp = reshape(res_X_north_cmp,grid_num_lon_pts, grid_num_lat_pts);
            res_Y_east_cmp = reshape(res_Y_east_cmp,grid_num_lon_pts, grid_num_lat_pts);
            res_Z_vert_cmp = reshape(res_Z_vert_cmp,grid_num_lon_pts, grid_num_lat_pts);

            res_X = res_X_north_cmp';
            res_Y = res_Y_east_cmp';
            res_Z = res_Z_vert_cmp';


            Bx = res_X;
            By = res_Y;
            Bz = res_Z;

            C_ALT = grid_height;
            C_LON = reshape(BX_VALS(:, 1),grid_num_lon_pts, grid_num_lat_pts);
            C_LON = C_LON';

            C_LAT = reshape(BX_VALS(:, 2),grid_num_lon_pts, grid_num_lat_pts);
            C_LAT = C_LAT';

            Bxx = Bx*0;
            Byx = Bx*0;
            Bzx = Bx*0;
            Bxy = Bx*0;
            Byy = Bx*0;
            Bzy = Bx*0;

            for i = 2 : grid_num_lat_pts -1
                for j = 2 : grid_num_lon_pts -1
                    vect_w = aux_vect_rot([Bx(i, j-1); By(i, j-1); Bz(i, j-1)], C_LON(i, j-1), C_LAT(i, j-1), C_LON(i, j), C_LAT(i, j));
                    vect_e = aux_vect_rot([Bx(i, j+1); By(i, j+1); Bz(i, j+1)], C_LON(i, j+1), C_LAT(i, j+1), C_LON(i, j), C_LAT(i, j));
                    vect_s = aux_vect_rot([Bx(i-1, j); By(i-1, j); Bz(i-1, j)], C_LON(i-1, j), C_LAT(i-1, j), C_LON(i, j), C_LAT(i, j));
                    vect_n = aux_vect_rot([Bx(i+1, j); By(i+1, j); Bz(i+1, j)], C_LON(i+1, j), C_LAT(i+1, j), C_LON(i, j), C_LAT(i, j));


                    vect_c = [Bx(i, j); By(i, j); Bz(i, j)];


                    r_w = C_ALT + EARTH_RADIUS_TESS_M;
                    lat_w = C_LAT(i, j-1)*pi/180.0;
                    lon_w = C_LON(i, j-1)*pi/180.0;

                    r_e = C_ALT + EARTH_RADIUS_TESS_M;
                    lat_e = C_LAT(i, j+1)*pi/180.0;
                    lon_e = C_LON(i, j+1)*pi/180.0;

                    r_s = C_ALT + EARTH_RADIUS_TESS_M;
                    lat_s = C_LAT(i-1, j)*pi/180.0;
                    lon_s = C_LON(i-1, j)*pi/180.0;

                    r_n = C_ALT + EARTH_RADIUS_TESS_M;
                    lat_n = C_LAT(i+1, j)*pi/180.0;
                    lon_n = C_LON(i+1, j)*pi/180.0;

                    r_c = C_ALT + EARTH_RADIUS_TESS_M;
                    lat_c = C_LAT(i, j)*pi/180.0;
                    lon_c = C_LON(i, j)*pi/180.0;


                    [x_e, y_e, z_e] = sph2cart(lon_e,lat_e,r_e);
                    [x_n, y_n, z_n] = sph2cart(lon_n,lat_n,r_n);


                        [x_w, y_w, z_w] = sph2cart(lon_w,lat_w,r_w);
                        ew_dist = norm([x_w, y_w, z_w]-[x_e, y_e, z_e]);

                        [x_s, y_s, z_s] = sph2cart(lon_s,lat_s,r_s);
                        ns_dist = norm([x_s, y_s, z_s]-[x_n, y_n, z_n]);

                        Bxyz_y = -(vect_w - vect_e)*1000000/(ew_dist) ;
                        Bxyz_x = -(vect_s - vect_n)*1000000/(ns_dist) ;

                        Bxx(i, j) = Bxyz_x(1);
                        Byx(i, j) = Bxyz_x(2);
                        Bzx(i, j) = Bxyz_x(3);

                        Bxy(i, j) = Bxyz_y(1);
                        Byy(i, j) = Bxyz_y(2);
                        Bzy(i, j) = Bxyz_y(3);

                end
            end


            Bxx = Bxx(2:end-1, 2:end-1);
            Byx = Byx(2:end-1, 2:end-1);
            Bzx = Bzx(2:end-1, 2:end-1);
            Bxy = Bxy(2:end-1, 2:end-1);
            Byy = Byy(2:end-1, 2:end-1);
            Bzy = Bzy(2:end-1, 2:end-1);

            TEST_SUM_Bxx = TEST_SUM_Bxx + Bxx;
            TEST_SUM_Byx = TEST_SUM_Byx + Byx;
            TEST_SUM_Bzx = TEST_SUM_Bzx + Bzx;
            TEST_SUM_Bxy = TEST_SUM_Bxy + Bxy;
            TEST_SUM_Byy = TEST_SUM_Byy + Byy;
            TEST_SUM_Bzy = TEST_SUM_Bzy + Bzy;

            Bxx = reshape(Bxx, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);
            Byx = reshape(Byx, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);
            Bzx = reshape(Bzx, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);
            Bxy = reshape(Bxy, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);
            Byy = reshape(Byy, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);
            Bzy = reshape(Bzy, (grid_num_lat_pts-2)*(grid_num_lon_pts-2), 1);

            A(:, ni)=[Bxx; Byx; Bzx; Bxy; Byy; Bzy]; 

        end
        waitbar(1,h,sprintf(strcat(BAR_TEXT, '\n100%%')));
        close(h);

        fileID = fopen(strcat(PATH_OUTPUT_BLOCK, 'design_matrix_grad.bin'),'w');
        fwrite(fileID,A,'float32');
        fclose(fileID); 

    end
    
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 3 finished\n');
end


%% %%%%%%%%%%%%%%%%%   STEP 4   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if STEP4
    close all;
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('STEP 4\n');
    
    
    if(~(exist(OBSERVED_DATA_FILENAME_BX) && exist(OBSERVED_DATA_FILENAME_BY) && exist(OBSERVED_DATA_FILENAME_BZ)))
        fprintf('Observed data files not found. Check filenames. \n');
        return
    end
        
    

    %crease inversion result folder
    if(~exist(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT)))
        mkdir(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT));
    end
    
    if (RESULT_FOLDER_VECT(end) ~= '\')
        RESULT_FOLDER_VECT = strcat(RESULT_FOLDER_VECT, '\');
        fprintf('Missing \\ added to the end of RESULT_FOLDER_VECT. \n');
    end

    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);
    N_LAYERS = N_FILES;

    x_0 = [];
    lon_c = [];
    lat_c = [];
    N_BODIES = [];
    A = [];

    iter_i = 0;
    for layer_i = 1 : N_LAYERS
        current = FILENAMES(layer_i).name;
        DATASET = sprintf('%s', current(1:end-14));
        sprintf('Current folder: %s\n',DATASET);

        LOCAL_OBSERVED_NAME = CURRENT_DATE;

        PATH = strcat(PATH_OUTPUT, DATASET,'\');

        BODIES = dir(fullfile(strcat(PATH, '*.magtess')));
        N_BODIES = [N_BODIES, length(BODIES)];


        dummy = textread(strcat(PATH, BODIES(1).name, '_Bz.txt'),'%f', 'commentstyle', 'shell');
        N_POINTS = length(dummy)/4;


        SUSCEPT_GRID_AVAILABLE(layer_i) = 0;

        %here you ust specify 
        APRIORI_SUSCEPT_FILENAME = sprintf(APRIORI_SUSCEPT_FILENAME_TEMPLATE, layer_i);
        if exist(APRIORI_SUSCEPT_FILENAME, 'file')
            [LON_s, LAT_s, SUSCEPT] = surface_read_from_file(APRIORI_SUSCEPT_FILENAME);
            F_s = scatteredInterpolant(LON_s,LAT_s, SUSCEPT, 'linear', 'none');
            SUSCEPT_GRID_AVAILABLE(layer_i) = 1;
        else
            F_s = NaN;
            fprintf('No a-priori suscept file %s\n', APRIORI_SUSCEPT_FILENAME);
        end

        suscept_init_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_layer', num2str(layer_i),'_suscept_init.xyz'), 'w');
        lon_c = [lon_c; zeros(N_BODIES(layer_i), 1)];
        lat_c = [lat_c; zeros(N_BODIES(layer_i), 1)];
        x_0 = [x_0; zeros(N_BODIES(layer_i), 1)];
        fprintf('Total tesseroids: %d\n', N_BODIES(layer_i));

        fileID = fopen(strcat(PATH, 'design_matrix.bin'),'r');
        A = [A, fread(fileID,[N_POINTS*3, N_BODIES(layer_i)],'float32')];
        fclose(fileID);

        BAR_TEXT = 'Loading susceptibilities...';
        h = waitbar(0,BAR_TEXT);
        for i = 1 : N_BODIES(layer_i)
                %PROGRESS BAR
            if rem(i, 100) == 0
                perc = (i/N_BODIES(layer_i));
                waitbar(perc,h,sprintf(strcat(BAR_TEXT, '\n%.1f%%'),perc*100));
            end

            one_tess = textread(strcat(PATH, BODIES(i).name),'%f', 'commentstyle', 'shell');
            lon_c(iter_i+i) = (one_tess(1) + one_tess(2)) / 2.0;
            lat_c(iter_i+i) = (one_tess(3) + one_tess(4)) / 2.0;
            if SUSCEPT_GRID_AVAILABLE(layer_i) == 0
                x_0(iter_i+i) = NaN;
            else
                x_0(iter_i+i) = F_s(lon_c(iter_i+i), lat_c(iter_i+i));
            end

            if (isnan(x_0(iter_i+i)))
                x_0(iter_i+i) = STANDART_SUSCEPTIBILITY(layer_i);
            end

            fprintf(suscept_init_fid, '%.2f %.2f %.12f\n', lon_c(iter_i+i), lat_c(iter_i+i), x_0(iter_i+i));     
        end
        waitbar(1,h,sprintf(strcat(BAR_TEXT, '\n100%%')));
        close(h);
        fclose(suscept_init_fid);

        iter_i = iter_i + N_BODIES(layer_i);
    end


    %% LOAD GRID
    current = FILENAMES(1).name;
    DATASET = sprintf('%s', current(1:end-14));
    sprintf('Current folder: %s\n',DATASET);


    PATH = strcat(PATH_OUTPUT, DATASET,'\');
    BODIES = dir(fullfile(strcat(PATH, '*.magtess')));

    BZ_GRID_FILENAME = strcat(PATH, BODIES(1).name, '_Bz.txt');
    BZ_GRID = textread(BZ_GRID_FILENAME,'%f', 'commentstyle', 'shell');
    N_POINTS = length(BZ_GRID)/4;
    BZ_GRID_VALS = reshape(BZ_GRID, 4, N_POINTS)';

    LON_GRID = BZ_GRID_VALS(:, 1);
    LAT_GRID = BZ_GRID_VALS(:, 2);


    %% LOAD OBSERVED DATA
    [LON, LAT, BX] = aux_surface_read_from_file(OBSERVED_DATA_FILENAME_BX);
    [LON, LAT, BY] = aux_surface_read_from_file(OBSERVED_DATA_FILENAME_BY);
    [LON, LAT, BZ] = aux_surface_read_from_file(OBSERVED_DATA_FILENAME_BZ);
    FX = scatteredInterpolant(LON, LAT, BX);
    FY = scatteredInterpolant(LON, LAT, BY);
    FZ = scatteredInterpolant(LON, LAT, BZ);

    d = zeros(N_POINTS*3, 1);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_observed_Bx.xyz'), 'w');
    for i = 1 : N_POINTS
        d(i) = FX(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_observed_By.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS+i) = FY(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_observed_Bz.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS*2+i) = -FZ(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), -d(N_POINTS*2+i));
    end
    fclose(datafid);

    original_d = d;


    %% INVERSION

    %% COVARIANCE MATRIX SUSCEPT

    COV_X = [];
    for di = 1 : N_LAYERS 
        RRR = diag(ones(N_BODIES(di), 1))*(sigma_x(di)^2);
        COV_X = blkdiag(COV_X, RRR);
    end

    Q = inv(COV_X);

    %% COVARIANCE MATRIX DATA

    COV_d = [];
    for di = 1 : 3
        COV_d = blkdiag(COV_d, diag(ones(N_POINTS, 1))*(sigma_vect_d(di)^2));
    end

    P = inv(COV_d);

    %% LEAST SQUARES

    fprintf('\nInversion:\n');

    L = (A.' * P * A + Q);
    R = (A.' * P * d + Q * x_0);
    chi = L \ R;

    %% SAVING RESULTING SUSCEPTIBILITIES TO FILES
    curr_i = 0;
    for layer_i = 1 : N_LAYERS
        suscapt_result_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_layer', num2str(layer_i),'_suscept_rslt.xyz'), 'w');

        for l = (curr_i+1) : (curr_i+N_BODIES(layer_i))
            fprintf(suscapt_result_fid, '%.2f %.2f %.12f\n', lon_c(l), lat_c(l), chi(l));  
        end

        fclose(suscapt_result_fid);
        curr_i = curr_i + N_BODIES(layer_i);
    end

    fprintf('\n');

    %% SAVING DIFFERENCE TO FILE S   
    result_field = aux_calc_anomaly_nograd_fast_allcomp(PATH_OUTPUT, DATASET, A, chi, strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME));

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_difference_Bx.xyz'), 'w');
    for i = 1 : N_POINTS
        dX = d(i)-result_field(i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dX);

    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_difference_By.xyz'), 'w');
    for i = 1 : N_POINTS
        dY = d(N_POINTS+i)-result_field(N_POINTS+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dY);
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_difference_Bz.xyz'), 'w');
    for i = 1 : N_POINTS
        dZ = -d(N_POINTS*2+i)+result_field(N_POINTS*2+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dZ);
    end
    fclose(datafid);


    %%make a picture
    WEST = WEST_EDGE;              
    EAST = EAST_EDGE;              
    SOUTH = SOUTH_EDGE;            
    NORTH = NORTH_EDGE; 
    
    SUSMIN = min(max(x_0), min(chi));
    SUSMAX = max(max(x_0), max(chi));
    STEP = round(((SUSMAX - SUSMIN)/10)*100)/100;

    fid_bat=fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, 'exec.bat'), 'w');

    for lay = 1 : N_LAYERS
        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_layer', num2str(lay),'_suscept_rslt'));
        fprintf(fid_bat, 'nearneighbor  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.5d -S1.5d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Cpanoply -E19 -T= -L%f/%f>  %%grid1%%.cpt\n', SUSMIN, SUSMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"Susceptibility":/:SI: -O >> %%grid1%%.ps\n', STEP);
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');


        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_layer', num2str(lay),'_suscept_init'));
        fprintf(fid_bat, 'nearneighbor  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.5d -S1.5d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Cpanoply -E19 -T= -L%f/%f >  %%grid1%%.cpt\n', SUSMIN, SUSMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"Susceptibility":/:SI: -O >> %%grid1%%.ps\n', STEP);
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');
    end

    FIELDMIN = min([min(-original_d), min(-d)]);
    FIELDMAX = max([max(-original_d), max(-d)]);
    FIELDSTEP = round(round(((FIELDMAX - FIELDMIN)/10)*100)/100);

    COMP_let = cellstr(['Bx'; 'By'; 'Bz']);
    for cl = 1 : length(COMP_let)

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, strcat('_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30 -L%f/%f >  %%grid1%%.cpt\n', FIELDMIN, FIELDMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"forward after inv. %s":/:nT: -O >> %%grid1%%.ps\n', FIELDSTEP, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, strcat('_observed_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30 -L%f/%f  >  %%grid1%%.cpt\n', FIELDMIN, FIELDMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"original observed %s":/:nT: -O >> %%grid1%%.ps\n', FIELDSTEP, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, strcat('_difference_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30  >  %%grid1%%.cpt\n');
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"(or. obs- - surr.)-(forward after inv.)%s":/:nT: -O >> %%grid1%%.ps\n', round(min(sigma_vect_d)*100)/100, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

    end  

    fclose(fid_bat);

    system(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, 'exec.bat'));

    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, '*.ps'));
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, '*.grd'));
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, '*.cpt'))
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, '*.d'))


    %% SAVING RESULTS

    FILNAME_MASK = strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, CURRENT_DATE);

    X0_FILENAME = strcat(FILNAME_MASK,'_x_0.bin');
    fileID = fopen(X0_FILENAME,'w');
    fwrite(fileID,x_0,'float32');
    fclose(fileID); 

    CHI0_FILENAME = strcat(FILNAME_MASK,'_chi.bin');
    fileID = fopen(CHI0_FILENAME,'w');
    fwrite(fileID,chi,'float32');
    fclose(fileID); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);

    fprintf('Files to process:\n');
    for n = 1 : N_FILES
        fprintf('  %d. %s\n', n, FILENAMES(n).name);
    end

    EARTH_RADIUS_TESS_M = 6378137.0;
    ALTITUDE = 400000; 

    fprintf('Calculating...\n')

    result_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_VECT, LOCAL_OBSERVED_NAME, '_RESULT_model.magtess'), 'w');

    for i = 1 : N_FILES
        current = FILENAMES(i).name;
        fid = fopen(strcat(PATH_OUTPUT, current), 'r');

        j = 1;
        tline = fgetl(fid);
        while ischar(tline)
            numbers = str2num(tline);

            fprintf(result_fid, '%d %d %d %d %f %f %f %f %.20f %.20f %.20f\n', numbers(1), numbers(2), numbers(3), numbers(4), numbers(5), numbers(6), numbers(7), chi(j), numbers(9), numbers(10), numbers(11));
            j = j + 1;


            tline = fgetl(fid);    
        end

        fclose(fid);


    end

    fclose(result_fid);
    
    
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 4 finished\n');
end


%% %%%%%%%%%%%%%%%%%   STEP 5   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if STEP5
    close all;
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('STEP 5\n');
    
    if(~(exist(OBSERVED_DATA_FILENAME_BX) && exist(OBSERVED_DATA_FILENAME_BY) && exist(OBSERVED_DATA_FILENAME_BZ)))
        fprintf('Observed data files not found. Check filenames. \n');
        return
    end
    
    %crease inversion result folder
    if(~exist(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD)))
        mkdir(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD));
    end
    
    if (RESULT_FOLDER_GRAD(end) ~= '\')
        RESULT_FOLDER_GRAD = strcat(RESULT_FOLDER_GRAD, '\');
        fprintf('Missing \\ added to the end of RESULT_FOLDER_GRAD. \n');
    end


    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);
    N_LAYERS = N_FILES;

    x_0 = [];
    lon_c = [];
    lat_c = [];
    N_BODIES = [];
    A = [];

    iter_i = 0;
    for layer_i = 1 : N_LAYERS
        current = FILENAMES(layer_i).name;
        DATASET = sprintf('%s', current(1:end-14));
        sprintf('Current folder: %s\n',DATASET);

        WEST = WEST_EDGE;              
        EAST = EAST_EDGE;              
        SOUTH = SOUTH_EDGE;            
        NORTH = NORTH_EDGE; 

        LOCAL_OBSERVED_NAME = CURRENT_DATE;
         
        PATH = strcat(PATH_OUTPUT, DATASET,'\');

        BODIES = dir(fullfile(strcat(PATH, '*.magtess')));
        N_BODIES = [N_BODIES, length(BODIES)];

        dummy = textread(strcat(PATH, BODIES(1).name, '_Bz.txt'),'%f', 'commentstyle', 'shell');
        grid_num_pts = length(dummy)/4;

        dummy_VALS = reshape(dummy, 4, grid_num_pts)';

        grid_west_border = min(dummy_VALS(:, 1));
        grid_east_border = max(dummy_VALS(:, 1));

        for grid_num_lon_pts = 2 : grid_num_pts
            if(dummy_VALS(1, 1) == dummy_VALS(grid_num_lon_pts, 1))
                grid_num_lon_pts = grid_num_lon_pts - 1;
                break;
            end
        end

        grid_num_lat_pts = grid_num_pts / grid_num_lon_pts;

        N_POINTS = (grid_num_lat_pts-2)*(grid_num_lon_pts-2)% * 6

        SUSCEPT_GRID_AVAILABLE(layer_i) = 0;
        APRIORI_SUSCEPT_FILENAME = sprintf(APRIORI_SUSCEPT_FILENAME_TEMPLATE, layer_i);
        if exist(APRIORI_SUSCEPT_FILENAME, 'file')
            [LON_s, LAT_s, SUSCEPT] = surface_read_from_file(APRIORI_SUSCEPT_FILENAME);
            F_s = scatteredInterpolant(LON_s,LAT_s, SUSCEPT, 'linear', 'none');
            SUSCEPT_GRID_AVAILABLE(layer_i) = 1;
        else
            F_s = NaN;
            fprintf('No apriori suscept file %s\n', APRIORI_SUSCEPT_FILENAME);
        end

        suscept_init_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_layer', num2str(layer_i),'_suscept_init.xyz'), 'w');
        lon_c = [lon_c; zeros(N_BODIES(layer_i), 1)];
        lat_c = [lat_c; zeros(N_BODIES(layer_i), 1)];
        x_0 = [x_0; zeros(N_BODIES(layer_i), 1)];
        fprintf('Total tesseroids: %d\n', N_BODIES(layer_i));

        fileID = fopen(strcat(PATH, 'design_matrix_grad.bin'),'r');
        A = [A, fread(fileID,[N_POINTS*6, N_BODIES(layer_i)],'float32')];
        fclose(fileID);  

        BAR_TEXT = 'Loading susceptibilities...';
        h = waitbar(0,BAR_TEXT);
        for i = 1 : N_BODIES(layer_i)
                %PROGRESS BAR
            if rem(i, 100) == 0
                perc = (i/N_BODIES(layer_i));
                waitbar(perc,h,sprintf(strcat(BAR_TEXT, '\n%.1f%%'),perc*100));
            end

            one_tess = textread(strcat(PATH, BODIES(i).name),'%f', 'commentstyle', 'shell');
            lon_c(iter_i+i) = (one_tess(1) + one_tess(2)) / 2.0;
            lat_c(iter_i+i) = (one_tess(3) + one_tess(4)) / 2.0;
            if SUSCEPT_GRID_AVAILABLE(layer_i) == 0
                x_0(iter_i+i) = NaN;
            else
                x_0(iter_i+i) = F_s(lon_c(iter_i+i), lat_c(iter_i+i));
            end

            if (isnan(x_0(iter_i+i)))
                x_0(iter_i+i) = STANDART_SUSCEPTIBILITY(layer_i);
            end

            fprintf(suscept_init_fid, '%.2f %.2f %.12f\n', lon_c(iter_i+i), lat_c(iter_i+i), x_0(iter_i+i));     
        end
        waitbar(1,h,sprintf(strcat(BAR_TEXT, '\n100%%')));
        close(h);
        fclose(suscept_init_fid);

        iter_i = iter_i + N_BODIES(layer_i);
    end

    %% LOAD GRID
    %OBSERVED DATA

    OBSERVED_BX = textread(OBSERVED_DATA_FILENAME_BX,'%f', 'commentstyle', 'shell');
    OBSERVED_BY = textread(OBSERVED_DATA_FILENAME_BY,'%f', 'commentstyle', 'shell');
    OBSERVED_BZ = textread(OBSERVED_DATA_FILENAME_BZ,'%f', 'commentstyle', 'shell');
    obs_grid_num_pts = length(OBSERVED_BX)/3;

    dummy_VALS = reshape(OBSERVED_BX , 3, obs_grid_num_pts)';

    obs_grid_west_border = min(dummy_VALS(:, 1));
    obs_grid_east_border = max(dummy_VALS(:, 1));

    for obs_grid_num_lon_pts = 2 : obs_grid_num_pts
        if(dummy_VALS(1, 1) ~= dummy_VALS(obs_grid_num_lon_pts, 1))
            obs_grid_num_lon_pts = obs_grid_num_lon_pts - 1;
            break;
        end
    end

    obs_grid_lat_step = abs(dummy_VALS(1, 2) - dummy_VALS(2, 2));

    obs_grid_num_lat_pts = obs_grid_num_pts / obs_grid_num_lon_pts;

    obs_grid_lon_step = abs(dummy_VALS(1, 1) - dummy_VALS(1+obs_grid_num_lat_pts, 1));

    BX_VALS = reshape(OBSERVED_BX, 3, obs_grid_num_pts)';
    BY_VALS = reshape(OBSERVED_BY, 3, obs_grid_num_pts)';
    BZ_VALS = reshape(OBSERVED_BZ, 3, obs_grid_num_pts)';

    obs_Bx = reshape(BX_VALS(:, 3),obs_grid_num_lon_pts, obs_grid_num_lat_pts);%
    obs_By = reshape(BY_VALS(:, 3),obs_grid_num_lon_pts, obs_grid_num_lat_pts);%
    obs_Bz = reshape(BZ_VALS(:, 3),obs_grid_num_lon_pts, obs_grid_num_lat_pts);%

    C_ALT = ALTITUDE;

    C_LON = reshape(BX_VALS(:, 1),obs_grid_num_lon_pts, obs_grid_num_lat_pts);
    C_LON = C_LON;%

    C_LAT = reshape(BX_VALS(:, 2),obs_grid_num_lon_pts, obs_grid_num_lat_pts);
    C_LAT = C_LAT;%

    obs_Bxx = obs_Bx*0;
    obs_Byx = obs_Bx*0;
    obs_Bzx = obs_Bx*0;
    obs_Bxy = obs_Bx*0;
    obs_Byy = obs_Bx*0;
    obs_Bzy = obs_Bx*0;

    for i = 2 : obs_grid_num_lon_pts -1
         for j = 2 : obs_grid_num_lat_pts -1
               vect_w = aux_vect_rot([obs_Bx(i, j-1); obs_By(i, j-1); obs_Bz(i, j-1)], C_LON(i, j-1), C_LAT(i, j-1), C_LON(i, j), C_LAT(i, j));
               vect_e = aux_vect_rot([obs_Bx(i, j+1); obs_By(i, j+1); obs_Bz(i, j+1)], C_LON(i, j+1), C_LAT(i, j+1), C_LON(i, j), C_LAT(i, j));
               vect_s = aux_vect_rot([obs_Bx(i-1, j); obs_By(i-1, j); obs_Bz(i-1, j)], C_LON(i-1, j), C_LAT(i-1, j), C_LON(i, j), C_LAT(i, j));
               vect_n = aux_vect_rot([obs_Bx(i+1, j); obs_By(i+1, j); obs_Bz(i+1, j)], C_LON(i+1, j), C_LAT(i+1, j), C_LON(i, j), C_LAT(i, j));

               vect_c = [obs_Bx(i, j); obs_By(i, j); obs_Bz(i, j)];

               r_w = C_ALT + EARTH_RADIUS_TESS_M;
               lat_w = C_LAT(i, j-1)*pi/180.0;
               lon_w = C_LON(i, j-1)*pi/180.0;

               r_e = C_ALT + EARTH_RADIUS_TESS_M;
               lat_e = C_LAT(i, j+1)*pi/180.0;
               lon_e = C_LON(i, j+1)*pi/180.0;

               r_s = C_ALT + EARTH_RADIUS_TESS_M;
               lat_s = C_LAT(i-1, j)*pi/180.0;
               lon_s = C_LON(i-1, j)*pi/180.0;

               r_n = C_ALT + EARTH_RADIUS_TESS_M;
               lat_n = C_LAT(i+1, j)*pi/180.0;
               lon_n = C_LON(i+1, j)*pi/180.0;

               r_c = C_ALT + EARTH_RADIUS_TESS_M;
               %lat_ic = C_LAT(i, j)*pi/180.0;
               %lon_ic = C_LON(i, j)*pi/180.0;


               [x_e, y_e, z_e] = sph2cart(lon_e,lat_e,r_e);
               [x_n, y_n, z_n] = sph2cart(lon_n,lat_n,r_n);

               [x_w, y_w, z_w] = sph2cart(lon_w,lat_w,r_w);
               ew_dist = norm([x_w, y_w, z_w]-[x_e, y_e, z_e]);

               [x_s, y_s, z_s] = sph2cart(lon_s,lat_s,r_s);
               ns_dist = norm([x_s, y_s, z_s]-[x_n, y_n, z_n]);

               Bxyz_y = -(vect_w - vect_e)*1000000/(ew_dist) ;
               Bxyz_x = -(vect_s - vect_n)*1000000/(ns_dist) ;

               obs_Bxx(i, j) = Bxyz_x(1);
               obs_Byx(i, j) = Bxyz_x(2);
               obs_Bzx(i, j) = Bxyz_x(3);

               obs_Bxy(i, j) = Bxyz_y(1);
               obs_Byy(i, j) = Bxyz_y(2);
               obs_Bzy(i, j) = Bxyz_y(3);

          end
    end

    obs_Bxx = obs_Bxx(2:end-1, 2:end-1);
    obs_Byx = obs_Byx(2:end-1, 2:end-1);
    obs_Bzx = obs_Bzx(2:end-1, 2:end-1);
    obs_Bxy = obs_Bxy(2:end-1, 2:end-1);
    obs_Byy = obs_Byy(2:end-1, 2:end-1);
    obs_Bzy = obs_Bzy(2:end-1, 2:end-1);

%{
            figure('Name', strcat('Gradients (calculated using central differences)'))
            sp_bxx = subplot(3, 3, 1);
            pcolor(obs_Bxx);
            shading flat
            shading interp
            cb_bxx = colorbar;
            xlabel(cb_bxx, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_bxx,'B_{xx}')

            sp_byx = subplot(3, 3, 4);
            pcolor(obs_Byx);
            shading flat
            shading interp
            cb_byx = colorbar;
            xlabel(cb_byx, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_byx,'B_{yx}')

            sp_bzx = subplot(3, 3, 7);
            pcolor(obs_Bzx);
            shading flat
            shading interp
            cb_bzx = colorbar;
            xlabel(cb_bzx, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_bzx,'B_{zx}')

            sp_bxy = subplot(3, 3, 2);
            pcolor(obs_Bxy);
            shading flat
            shading interp
            cb_bxy = colorbar;
            xlabel(cb_bxy, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_bxy,'B_{xy}');

            sp_byy = subplot(3, 3, 5);
            pcolor(obs_Byy);
            shading flat
            shading interp
            cb_byy = colorbar;
            xlabel(cb_byy, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_byy,'B_{yy}');


            sp_bzy = subplot(3, 3, 8);
            pcolor(obs_Bzy);
            shading flat
            shading interp
            cb_bzy = colorbar;
            xlabel(cb_bzy, '[pT/km]');
            axis tight;
            axis square;
            ylabel('latitude [deg]');
            xlabel('longitude [deg]');
            zlabel('[pT]/km');
            title(sp_bzy,'B_{zy}');
%}

    C_LON = C_LON(2:end-1, 2:end-1);
    C_LAT = C_LAT(2:end-1, 2:end-1);


    obs_Bxx = reshape(obs_Bxx, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Byx = reshape(obs_Byx, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Bzx = reshape(obs_Bzx, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Bxy = reshape(obs_Bxy, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Byy = reshape(obs_Byy, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Bzy = reshape(obs_Bzy, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    obs_Bzz = -obs_Bxx -obs_Byy;

    C_LON = reshape(C_LON, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);
    C_LAT = reshape(C_LAT, (obs_grid_num_lat_pts-2)*(obs_grid_num_lon_pts-2), 1);

    FXX = scatteredInterpolant(C_LON, C_LAT, obs_Bxx);
    FYX = scatteredInterpolant(C_LON, C_LAT, obs_Byx);
    FZX = scatteredInterpolant(C_LON, C_LAT, obs_Bzx);
    FXY = scatteredInterpolant(C_LON, C_LAT, obs_Bxy);
    FYY = scatteredInterpolant(C_LON, C_LAT, obs_Byy);
    FZY = scatteredInterpolant(C_LON, C_LAT, obs_Bzy);
    FZZ = scatteredInterpolant(C_LON, C_LAT, obs_Bzz);


    %%%%%%%%%%%%%%%
    current = FILENAMES(1).name;
    DATASET = sprintf('%s', current(1:end-14));

    PATH = strcat(PATH_OUTPUT, DATASET,'\');
    BODIES = dir(fullfile(strcat(PATH, '*.magtess')));

    BZ_GRID_FILENAME = strcat(PATH, BODIES(1).name, '_Bz.txt');
    BZ_GRID = textread(BZ_GRID_FILENAME,'%f', 'commentstyle', 'shell');
    tN_POINTS = length(BZ_GRID)/4;
    BZ_GRID_VALS = reshape(BZ_GRID, 4, tN_POINTS)';

    LON_GRID = BZ_GRID_VALS(:, 1);
    LAT_GRID = BZ_GRID_VALS(:, 2);

    LON_GRID = reshape(LON_GRID,grid_num_lon_pts, grid_num_lat_pts)';
    LAT_GRID = reshape(LAT_GRID,grid_num_lon_pts, grid_num_lat_pts)';


    LON_GRID = LON_GRID(2:end-1, 2:end-1);
    LAT_GRID = LAT_GRID(2:end-1, 2:end-1);


    LON_GRID = reshape(LON_GRID,(grid_num_lon_pts -2) * (grid_num_lat_pts-2), 1);
    LAT_GRID = reshape(LAT_GRID,(grid_num_lon_pts -2) * (grid_num_lat_pts-2), 1);



    d = zeros(N_POINTS*6, 1);%%%

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Bxx.xyz'), 'w');
    for i = 1 : N_POINTS
        d(i) = FXX(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Byx.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS+i) = FYX(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Bzx.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS*2+i) = FZX(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS*2+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Bxy.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS*3+i) = FXY(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS*3+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Byy.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS*4+i) = FYY(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS*4+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_observed_Bzy.xyz'), 'w');
    for i = 1 : N_POINTS
        d(N_POINTS*5+i) = FZY(LON_GRID(i), LAT_GRID(i));
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), d(N_POINTS*5+i));
    end
    fclose(datafid);

    original_d = d;


    %% COVARIANCE MATRIX SUSCEPT

    COV_X = [];
    for di = 1 : N_LAYERS
        RRR = diag(ones(N_BODIES(di), 1))*(sigma_x(di)^2);

        COV_X = blkdiag(COV_X, RRR);
    end

    Q = inv(COV_X);

    %% COVARIANCE MATRIX DATA

    COV_d = [];
    for di = 1 : 6
        COV_d = blkdiag(COV_d, diag(ones(N_POINTS, 1))*(sigma_grad_d(di)^2));
    end

    P = inv(COV_d);


    fprintf('\nInversion:\n');


    L = (A.' * P * A + Q);
    R = (A.' * P * d + Q * x_0);
    chi = L \ R;


    curr_i = 0;
    for layer_i = 1 : N_LAYERS
        suscapt_result_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_layer', num2str(layer_i),'_suscept_rslt.xyz'), 'w');

        for l = (curr_i+1) : (curr_i+N_BODIES(layer_i))
            fprintf(suscapt_result_fid, '%.2f %.2f %.12f\n', lon_c(l), lat_c(l), chi(l));  
        end

        fclose(suscapt_result_fid);
        curr_i = curr_i + N_BODIES(layer_i);
    end

    fprintf('\n');


    result_field = aux_calc_anomaly_nograd_fast_allcomp(PATH_OUTPUT, DATASET, A, chi, strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME));

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Bxx.xyz'), 'w');
    for i = 1 : N_POINTS
        dXX = d(i)-result_field(i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dXX);

    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Byx.xyz'), 'w');
    for i = 1 : N_POINTS
        dYX = d(N_POINTS+i)-result_field(N_POINTS+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dYX);
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Bzx.xyz'), 'w');
    for i = 1 : N_POINTS
        dZX = d(N_POINTS*2+i)-result_field(N_POINTS*2+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dZX);
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Bxy.xyz'), 'w');
    for i = 1 : N_POINTS
        dXY = d(N_POINTS*3+i)-result_field(N_POINTS*3+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dXY);

    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Byy.xyz'), 'w');
    for i = 1 : N_POINTS
        dYY = d(N_POINTS*4+i)-result_field(N_POINTS*4+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dYY);
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_difference_Bzy.xyz'), 'w');
    for i = 1 : N_POINTS
        dZY = d(N_POINTS*5+i)-result_field(N_POINTS*5+i);
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), dZY);
    end
    fclose(datafid);


    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Bxx.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i),  result_field(i));

    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Byx.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), result_field(N_POINTS+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Bzx.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), result_field(N_POINTS*2+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Bxy.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), result_field(N_POINTS*3+i));

    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Byy.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), result_field(N_POINTS*4+i));
    end
    fclose(datafid);

    datafid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_result_Bzy.xyz'), 'w');
    for i = 1 : N_POINTS
        fprintf(datafid, '%.2f %.2f %.12f\n', LON_GRID(i), LAT_GRID(i), result_field(N_POINTS*5+i));
    end
    fclose(datafid);

    %%make a picture
    WEST = WEST_EDGE;              
    EAST = EAST_EDGE;              
    SOUTH = SOUTH_EDGE;            
    NORTH = NORTH_EDGE; 

    SUSMIN = min(min(x_0), min(chi));
    SUSMAX = max(max(x_0), max(chi));
    STEP = round(((SUSMAX - SUSMIN)/10)*1000)/1000;



    fid_bat=fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, 'exec.bat'), 'w');

    for lay = 1 : N_LAYERS
        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_layer', num2str(lay),'_suscept_rslt'));
        fprintf(fid_bat, 'nearneighbor  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.5d -S1.5d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Cpanoply -E19 -T= -L%f/%f>  %%grid1%%.cpt\n', SUSMIN, SUSMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"Susceptibility":/:SI: -O >> %%grid1%%.ps\n', STEP);
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');


        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, '_layer', num2str(lay),'_suscept_init'));
        fprintf(fid_bat, 'nearneighbor  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.5d -S1.5d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Cpanoply -E19 -T= -L%f/%f >  %%grid1%%.cpt\n', SUSMIN, SUSMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psxy %s_edges.xy -Rd%d/%d/%d/%d -JU33N/6i -O -K -Wthickest,white,- >> %%grid1%%.ps\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME), WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"Susceptibility":/:SI: -O >> %%grid1%%.ps\n', STEP);
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');
    end

    FIELDMIN = min([min(-original_d), min(-d)]);
    FIELDMAX = max([max(-original_d), max(-d)]);
    FIELDSTEP = round(round(((FIELDMAX - FIELDMIN)/10)*100)/100);

    COMP_let = cellstr(['Bxx'; 'Byx'; 'Bzx'; 'Bxy'; 'Byy'; 'Bzy']);
    for cl = 1 : length(COMP_let)

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, strcat('_result_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30 -L%f/%f >  %%grid1%%.cpt\n', FIELDMIN, FIELDMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"forward after inv. %s":/:"pT/km": -O >> %%grid1%%.ps\n', FIELDSTEP, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, strcat('_observed_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30 -L%f/%f  >  %%grid1%%.cpt\n', FIELDMIN, FIELDMAX);
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"original observed %s":/:"pT/km": -O >> %%grid1%%.ps\n', FIELDSTEP, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

        fprintf(fid_bat, 'set grid1=%s\n', strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, LOCAL_OBSERVED_NAME, strcat('_difference_', char(COMP_let(cl))) ) );
        fprintf(fid_bat, 'surface  %%grid1%%.xyz -Rd%d/%d/%d/%d -I0.25d -G%%grid1%%.grd\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'grd2cpt %%grid1%%.grd -Chaxby -E30  >  %%grid1%%.cpt\n');
        fprintf(fid_bat, 'grdimage %%grid1%%.grd -Rd%d/%d/%d/%d  -JU33N/6i -C%%grid1%%.cpt -P -K > %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'pscoast -Dl -Rd%d/%d/%d/%d  -JU33N/6i -Ba5g5f1/a5g5f1WESN   -V -W1/0.5thin   -O -K >> %%grid1%%.ps\n', WEST, EAST, SOUTH, NORTH);
        fprintf(fid_bat, 'psscale -D5.3i/5i/5c/0.7c -C%%grid1%%.cpt -I -B%f:"(or. obs- - surr.)-(forward after inv.)%s":/:"pT/km": -O >> %%grid1%%.ps\n', round(min(sigma_grad_d)*100)/100, char(COMP_let(cl)));
        fprintf(fid_bat, 'ps2raster -A+r %%grid1%%.ps\n\n');

    end  

    fclose(fid_bat);

    system(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, 'exec.bat'));

    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, '*.ps'));
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, '*.grd'));
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, '*.cpt'))
    delete(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, '*.d'))
    %%%%%%%%%%%%%%%%%%%%

    FILNAME_MASK = strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, CURRENT_DATE);

    X0_FILENAME = strcat(FILNAME_MASK,'_x_0.bin');
    fileID = fopen(X0_FILENAME,'w');
    fwrite(fileID,x_0,'float32');
    fclose(fileID); 

    CHI0_FILENAME = strcat(FILNAME_MASK,'_chi.bin');
    fileID = fopen(CHI0_FILENAME,'w');
    fwrite(fileID,chi,'float32');
    fclose(fileID); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FILENAMES = dir(fullfile(strcat(PATH_OUTPUT, '*.magtess_block')));
    N_FILES = length(FILENAMES);

    fprintf('Files to process:\n');
    for n = 1 : N_FILES
        fprintf('  %d. %s\n', n, FILENAMES(n).name);
    end

    EARTH_RADIUS_TESS_M = 6378137.0;
    ALTITUDE = 400000; 

    fprintf('Calculating...\n')

    result_fid = fopen(strcat(PATH_OUTPUT, RESULT_FOLDER_GRAD, '_RESULT_model.magtess'), 'w');

    for i = 1 : N_FILES
        current = FILENAMES(i).name;

        fid = fopen(strcat(PATH_OUTPUT, current), 'r');


        j = 1;
        tline = fgetl(fid);
        while ischar(tline)
            numbers = str2num(tline);

            fprintf(result_fid, '%d %d %d %d %f %f %f %f %.20f %.20f %.20f\n', numbers(1), numbers(2), numbers(3), numbers(4), numbers(5), numbers(6), numbers(7), chi(j), numbers(9), numbers(10), numbers(11));
            j = j + 1;


            tline = fgetl(fid);    
        end

        fclose(fid);


    end

    fclose(result_fid);
    
    
    clearvars -except EARTH_RADIUS_IGRF_KM EARTH_RADIUS_TESS_M  STEP1 STEP2 STEP3 STEP4 STEP5 IGRF_DATE PATH_SURF TOP_FILENAME BOT_FILENAME PATH_OUTPUT TESSEROID_WIDTH LAYER_HEIGHT WEST_EDGE EAST_EDGE SOUTH_EDGE NORTH_EDGE EDGE_EXTENSION_W EDGE_EXTENSION_E EDGE_EXTENSION_S EDGE_EXTENSION_N ALTITUDE SPACING OBSERVED_DATA_FILENAME_BX OBSERVED_DATA_FILENAME_BY OBSERVED_DATA_FILENAME_BZ RESULT_FOLDER_VECT RESULT_FOLDER_GRAD STANDART_SUSCEPTIBILITY APRIORI_SUSCEPT_FILENAME_TEMPLATE sigma_x sigma_vect_d sigma_grad_d BLOCK_FILENAME_TEMPLATE CURRENT_DATE
    fprintf('STEP 5 finished\n');
end