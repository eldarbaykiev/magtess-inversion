%function D_tikhonov_for_block()
clc;
clear all;


MATLAB_FOLDER = '\\Filer2\Baykiev_Eldar\Dokumenter\MATLAB\';
cd(MATLAB_FOLDER);

PATH_DATA = 'C:\INVERSION_ALL_LAYERS\'
EDGE = 5;


FILENAMES = dir(fullfile(strcat(PATH_DATA, '*.magtess_block')));
N_FILES = length(FILENAMES);

N_BODIES = []
for layer_i = 1 : N_FILES
current = FILENAMES(layer_i ).name;
    fprintf('Block: %s\n', current);

    borders = sscanf(current, 'allblock_G1_W%dE%dS%dN%d');
    WEST = borders(1)+EDGE;              
    EAST = borders(2)-EDGE;              
    SOUTH = borders(3)+EDGE;            
    NORTH = borders(4)-EDGE;

    layer_depth = regExtractNums(current); %can be corrected
    HOT = -layer_depth(7);
    HOB = -layer_depth(6);
   
    DATASET = sprintf('W%dE%dS%dN%d_layer_%d_%d', WEST-EDGE, EAST+EDGE, SOUTH-EDGE, NORTH+EDGE, HOB, HOT);


    PATH = strcat(PATH_DATA, DATASET,'\');
        
    BODIES = dir(fullfile(strcat(PATH, '*.magtess')));
    N_BODIES = [N_BODIES, length(BODIES)];
    
    
end

layers = N_BODIES;
N_BODIES = sum(N_BODIES)

close all


PATH_TO_RESULT = 'C:\INVERSION_ALL_LAYERS\SIFM_TRUNCATION_DEGREE_13\'
filename1 = '090616_160051_chi.bin1';

fileID = fopen(strcat(PATH_TO_RESULT, filename1),'r');
x1 = fread(fileID,[1, N_BODIES],'float32');
fclose(fileID);


FILENAMES = dir(fullfile(strcat(PATH_TO_RESULT, '*_chi.bin3')));
FILENAMES0 = dir(fullfile(strcat(PATH_TO_RESULT, '*_x_0.bin3')));
N_FILES = length(FILENAMES);

layer_end = 0;
for l = 1 : 6
    
    
layer_start = layer_end+1;
layer_end = layer_start+layers(l)-1;;


figure('Name', strcat('Layer', num2str(l)))
stri = 'rgbkm'

for i = 1 : 1
   filename2 = FILENAMES(i).name 
   
    
    
    filename3 = FILENAMES0(i).name 
    
    
    fileID = fopen(strcat(PATH_TO_RESULT, filename3),'r');
x3 = fread(fileID,[1, N_BODIES],'float32');
fclose(fileID);



fileID = fopen(strcat(PATH_TO_RESULT, filename2),'r');
x2 = fread(fileID,[1, N_BODIES],'float32');
fclose(fileID);



scatter(x1(layer_start: layer_end), x2(layer_start: layer_end), '.', stri(i));
xlabel(strcat('Suscept [SI], a priori setup 0.015 SI'))
ylabel(strcat('Suscept [SI]'))
hold on

line ([0 0.1], [0 0.1])

line ([0 0], [-0.04 0.1])
line ([-0.04 0.1], [0 0])

x3(1)
stri(i)

end

end