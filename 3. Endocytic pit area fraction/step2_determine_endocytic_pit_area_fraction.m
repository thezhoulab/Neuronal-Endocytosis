%% This step is used to analyze area fraction and intensity for endocytic pits

clear all
close all
clc

parent_folder='Y:\20241208 D54D2 in b2KD and WT_APP\Analysis';
dirs=dir(parent_folder);
for m=1%:length(dirs)
    if dirs(m).isdir && ~strcmp(dirs(m).name, '.') && ~strcmp(dirs(m).name, '..')
        folder_path = fullfile(parent_folder, dirs(m).name);
    end
    cellmask_list = dir(fullfile(folder_path,'step 1_cell mask_MAP2','cellmask_*.tif'));
    puncta_list = dir(fullfile(folder_path,'D54D2_Cy3_*.tif')); % adjust name
    saving_path = fullfile(folder_path,'D54D2_Cy3'); % adjust name
    mkdir(saving_path);
    file_path = fullfile(folder_path,'results_D54D2_Cy3.txt'); % adjust name
    fileID = fopen(file_path,'w');
    fprintf(fileID,'neuron\tintensity\tcount\tarea_fraction\tpuncta_intensity\tmean_puncta_size\n');
    
    for i = 1:length(cellmask_list)
        img = imread(fullfile(folder_path,'step 1_cell mask_MAP2',cellmask_list(i).name));
        img = im2double(img);
        puncta = imread(fullfile(folder_path,puncta_list(i).name));
        original = puncta;
        intensity = mean(puncta(img > 0)); 
        puncta = imbinarize(puncta,"adaptive",'sensitivity',0.1);
        puncta = im2double(puncta);
        puncta = puncta .* img;
        se = strel('disk',2);
        puncta = imopen(puncta, se);
        cc = bwconncomp(puncta);
        min_size = 10;
        for z = 1:cc.NumObjects
            if numel(cc.PixelIdxList{z}) < min_size
                puncta(cc.PixelIdxList{z}) = 0;
            end
        end
        puncta = imfill(puncta, 'holes');
        puncta = puncta .* img;
        puncta_intensity = mean(original(puncta > 0));
        area_fraction = bwarea(puncta) / bwarea(img);
        cc = bwconncomp(puncta);
        count = cc.NumObjects;
        puncta_sizes = regionprops(cc, 'Area');
        puncta_sizes = [puncta_sizes.Area];
        mean_puncta_size = mean(puncta_sizes);
        fprintf(fileID,'%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n', ...
            i, intensity, count, area_fraction, puncta_intensity, mean_puncta_size);
        
        img = im2double(img);
        puncta = im2double(puncta);
        puncta = cat(3, puncta, puncta, puncta);
        alpha = 0.5 * ones(size(puncta, 1), size(puncta, 2));
        merged = img;
        for j = 1:size(img, 3)
            merged(:, :, j) = img(:, :, j) .* (1 - alpha) + puncta(:, :, j) .* alpha;
        end
        
        figure;
        subplot(1,2,1), imshow(imadjust(original)), title('original');
        subplot(1,2,2), imshow(merged), title('cell mask');
        
        if i < 10
            filename = fullfile(saving_path, sprintf('D54D2_Cy3_000%d.tif', i));
        else
            filename = fullfile(saving_path, sprintf('D54D2_Cy3_00%d.tif', i));
        end
        imwrite(merged, filename);
    end
    fclose(fileID);
end
