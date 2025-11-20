%% This file is used to calculate experimental and simulated Pearson’s correlation coefficients for surface-localized endocytic pits

clear all;
close all;
clc;

% Define folder paths
folder_path = 'Y:\20250927 DIV14 OE_EndoA2_583 b2_647\Analysis\diat_regions\2 Type Puncta';
saving_path = fullfile(folder_path, 'PCC');
if ~exist(saving_path, 'dir')
    mkdir(saving_path);
end

% Get the list of TIFF files in the folder
file_list = dir(fullfile(folder_path, '*.tif'));

XC = zeros(1, length(file_list));
XC_simulation_average = zeros(1, length(file_list));

% Loop through each file in the folder
for i = 1:length(file_list)
    % Read the multi-channel TIFF image using BioformatsImage
    img_path = fullfile(folder_path, file_list(i).name);
    try
        reader = BioformatsImage(img_path);
    catch ME
        fprintf('Failed to read image %s: %s\n', file_list(i).name, ME.message);
        continue; % Skip to the next file
    end
    
    % Check if the image has at least three channels
    if reader.sizeC < 3
        fprintf('Image %s does not have three channels. Skipping.\n', file_list(i).name);
        continue; % Skip this file and continue with the next one
    end
    
    % Read the first two channels
    channel1 = getPlane(reader, 1, 1, 1); % Read the first channel (background)
    channel2 = getPlane(reader, 1, 2, 1); % Read the second channel (clusters)
    
    % Make each image square by padding
    [height, width] = size(channel1);
    if height > width
        pad_size = (height - width) / 2;
        channel1 = padarray(channel1, [0, floor(pad_size)], 0, 'pre');
        channel1 = padarray(channel1, [0, ceil(pad_size)], 0, 'post');
        channel2 = padarray(channel2, [0, floor(pad_size)], 0, 'pre');
        channel2 = padarray(channel2, [0, ceil(pad_size)], 0, 'post');
    elseif width > height
        pad_size = (width - height) / 2;
        channel1 = padarray(channel1, [floor(pad_size), 0], 0, 'pre');
        channel1 = padarray(channel1, [ceil(pad_size), 0], 0, 'post');
        channel2 = padarray(channel2, [floor(pad_size), 0], 0, 'pre');
        channel2 = padarray(channel2, [ceil(pad_size), 0], 0, 'post');
    end
    
    % Normalize and enhance contrast for Channel 1 and Channel 2
    channel1_contrast = imadjust(im2double(channel1), stretchlim(im2double(channel1), [0.01, 0.99]), []);
    channel2_contrast = imadjust(im2double(channel2), stretchlim(im2double(channel2), [0.01, 0.99]), []);
    
    % Perform image segmentation to identify clusters in Channel 2
    threshold_level = graythresh(channel2_contrast);
    binary_image = imbinarize(channel2_contrast, threshold_level);
    se = strel('disk', 4); % Structuring element with a small radius for minimal erosion
    binary_image = imerode(binary_image, se); % Erode the binary image to shrink the clusters slightly
    
    % Label clusters
    [labeled_clusters, num_clusters] = bwlabel(binary_image);
    cluster_properties = regionprops(labeled_clusters, 'Area', 'Centroid');
    
    fprintf('Processing file: %s - Found %d clusters.\n', file_list(i).name, num_clusters);
    
    % Initialize an empty image for selected clusters' intensities
    selected_clusters = zeros(size(channel2), 'like', channel2);  % To store selected clusters with original intensities
    
    % Process each cluster that meets the size criteria % optional
    for j = 1:num_clusters
        % Compute the area, major axis length, and minor axis length
        area = cluster_properties(j).Area;
        major_axis_length = regionprops(labeled_clusters == j, 'MajorAxisLength');
        minor_axis_length = regionprops(labeled_clusters == j, 'MinorAxisLength');
        
        % Check if both properties were calculated correctly
        if ~isempty(major_axis_length) && ~isempty(minor_axis_length)
            major_axis_length = major_axis_length.MajorAxisLength;
            minor_axis_length = minor_axis_length.MinorAxisLength;
            
            % Calculate the ellipticity
            ellipticity = 1 - (minor_axis_length / major_axis_length);
        else
            fprintf('Cluster %d has missing properties for ellipticity calculation. Skipping.\n', j);
            continue; % Skip to the next cluster if any property is missing
        end
        
        % Only consider clusters that meet the size, ellipticity, and overlap criteria
        if area >= 50 && area <= 400 && ellipticity <= 1 % area bettwen 50 and 400
            fprintf('Cluster %d meets the size and ellipticity criteria. Area: %f\n', j, area);
            
            % Create a binary mask for the j-th cluster
            cluster_mask = (labeled_clusters == j);
            
            % Cast the logical cluster_mask to the same type as channel2
            cluster_mask = cast(cluster_mask, class(channel2));
            
            % Retain the original intensity of the selected clusters in Channel 2
            selected_clusters = selected_clusters + (cluster_mask .* channel2);
            
        end
    end
    
    % If no clusters meet the criteria, skip saving the image
    if all(selected_clusters(:) == 0)
        fprintf('No clusters met the criteria in file %s.\n', file_list(i).name);
        continue;
    end
    
    % Normalize the selected cluster intensities
    selected_clusters_norm = imadjust(selected_clusters, stretchlim(selected_clusters, [0.01, 0.99]), []);
    
    
    % Crop the image based on non-zero region of `imbw_1`
    % Get the non-zero x, y extreme values from imbw_1
    [rows, cols] = find(channel1_contrast > 0);
    if ~isempty(rows) && ~isempty(cols)
        row_min = min(rows);
        row_max = max(rows);
        col_min = min(cols);
        col_max = max(cols);
        
        % Crop both images based on the extreme values
        cropped_channel1 = channel1_contrast(row_min:row_max, col_min:col_max);
        cropped_selected_clusters = selected_clusters_norm(row_min:row_max, col_min:col_max);
        
        % Create the final combined image using the cropped regions
        final_combined_image = cat(3, uint8(cropped_channel1 * 255), uint8(cropped_selected_clusters * 255), zeros(size(cropped_channel1), 'uint8'));
    else
        fprintf('No non-zero pixels found in channel 1 for cropping.\n');
        continue; % Skip to the next image if there's no non-zero region
    end
    
    % Extract the file name without the extension
    [~, baseFileName, ~] = fileparts(file_list(i).name);
    
    % Save the combined image
    output_file_path = fullfile(saving_path, [baseFileName, '_clusters_overlay.tif']);
    imwrite(final_combined_image, output_file_path);
    figure
    imshow(final_combined_image);
    
    imbw_1=cropped_channel1;
    imbw_2=cropped_selected_clusters;
    imbw_1=double(imbw_1);
    imbw_2=double(imbw_2);
    XC(i) = corr2(imbw_1, imbw_2); % Calculate corr2()
    
    fprintf('Saved combined cluster image for %s.\n', baseFileName);
    
    % Simulation part: shift clusters randomly and calculate XC_simulation
    imbw_2_simulation = bwlabel(imbw_2);
    [max_row, max_col] = size(imbw_2_simulation);
    [index_row_1, index_col_1] = find(imbw_1 == 1);
    for simulation_count = 1:100
        for m = 1:max(imbw_2_simulation(:))
            [index_row, index_col] = find(imbw_2_simulation == m);
            if isempty(index_row)
                continue;
            end
            index_1 = find(imbw_2_simulation == m);
            imbw_2_simulation(index_1) = 0;
            for iii = 1:10000
                new_vector = 0;
                random_vector_x = randi([-min(index_row) + 1, max(index_row_1) - max(index_row)]);
                random_vector_y = randi([-min(index_col) + 1, max(index_col_1) - max(index_col)]);
                for mm = 1:length(index_row)
                    if imbw_2_simulation(min(max(index_row(mm) + random_vector_x, 1), max_row), ...
                            min(max(index_col(mm) + random_vector_y, 1), max_col)) > 0
                        new_vector = 1;
                    end
                end
                if new_vector == 0
                    break;
                end
            end
            if min(index_row + random_vector_x) < 1
                random_vector_x = random_vector_x + 1 - min(index_row + random_vector_x);
            end
            if min(index_col + random_vector_y) < 1
                random_vector_y = random_vector_y + 1 - min(index_col + random_vector_y);
            end
            for mm = 1:length(index_row)
                imbw_2_simulation(min(max(index_row(mm) + random_vector_x, 1), max_row), ...
                    min(max(index_col(mm) + random_vector_y, 1), max_col)) = m;
            end
        end
        XC_simulation(simulation_count) = corr2(imbw_1, imbw_2_simulation);
    end
    XC_simulation_average(i) = mean(XC_simulation);
end

XC = XC';
XC_simulation_average = XC_simulation_average';
