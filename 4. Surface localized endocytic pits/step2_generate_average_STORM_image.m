%% This file is used to generate dual-color STORM images by aligning individual STORM images to the centers of endocytic pits

clear all;
close all;
clc;

% Define folder paths
folder_path = 'Y:\20250927 DIV14 OE_EndoA2_583 b2_647\Analysis\diat_regions\2 Type Puncta';
saving_path = fullfile(folder_path, '2D average');
if ~exist(saving_path, 'dir')
    mkdir(saving_path);
end

% Get the list of TIFF files in the folder
file_list = dir(fullfile(folder_path, '*.tif'));

% Initialize a list to store all rotated cluster images
all_rotated_images = {};

% Dimming factor for making channels dimmer
dimming_factor = 0.5; % Adjust this value between 0 and 1 to control the level of dimming

% Loop through each file in the folder
for i = 1:length(file_list)
    % Read the multi-channel TIFF image using BioformatsImage
    img_path = fullfile(folder_path, file_list(i).name);
    try
        reader = BioformatsImage(img_path);
    catch ME
        fprintf('Failed to read image %s: %s\n', file_list(i).name, ME.message);
        continue; 
    end
    
    % Check if the image has at least three channels
    if reader.sizeC < 3
        fprintf('Image %s does not have three channels. Skipping.\n', file_list(i).name);
        continue; % Skip this file and continue with the next one
    end
    
    % Read the first two channels
    channel1 = getPlane(reader, 1, 1, 1); % Read the first channel
    channel2 = getPlane(reader, 1, 2, 1); % Read the second channel
    
    % Make each image square by padding
    [height, width] = size(channel1);
    max_dim = max(height, width);
    new_size = ceil(sqrt(2) * max_dim); % Calculate new size for diagonal coverage
    
    % Pad images to make them square and large enough to prevent clipping during rotation
    pad_height = (new_size - height) / 2;
    pad_width = (new_size - width) / 2;
    channel1 = padarray(channel1, [floor(pad_height), floor(pad_width)], 0, 'pre');
    channel1 = padarray(channel1, [ceil(pad_height), ceil(pad_width)], 0, 'post');
    channel2 = padarray(channel2, [floor(pad_height), floor(pad_width)], 0, 'pre');
    channel2 = padarray(channel2, [ceil(pad_height), ceil(pad_width)], 0, 'post');
    
    % Normalize and enhance contrast for Channel 1 and Channel 2
    channel1_contrast = imadjust(im2double(channel1), stretchlim(im2double(channel1), [0.01, 0.99]), []);
    channel2_contrast = imadjust(im2double(channel2), stretchlim(im2double(channel2), [0.01, 0.99]), []);
    
    % Apply dimming to both channels
    channel1_dimmer = channel1_contrast * dimming_factor;
    channel2_dimmer = channel2_contrast * dimming_factor;
    
    % Perform image segmentation to identify clusters in Channel 2
    threshold_level = graythresh(channel2_contrast);
    binary_image = imbinarize(channel2_contrast, threshold_level);
    se = strel('disk', 4); 
    binary_image = imerode(binary_image, se); 
    
    % Label clusters
    [labeled_clusters, num_clusters] = bwlabel(binary_image);
    cluster_properties = regionprops(labeled_clusters, 'Area', 'Centroid');
    
    fprintf('Processing file: %s - Found %d clusters.\n', file_list(i).name, num_clusters);
    
    % Process each cluster that meets the size criteria
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
            continue; 
        end
        
        if cluster_properties(j).Area >= 50 && cluster_properties(j).Area <= 400 && ellipticity <=1
            fprintf('Cluster %d meets the size criteria. Area: %f\n', j, cluster_properties(j).Area);
            
            % Create a binary mask for the j-th cluster
            cluster_mask = (labeled_clusters == j);
            
            % Cast the logical mask to the same class as channel2
            cluster_mask = cast(cluster_mask, class(channel2));
            
            % Retain the original intensity of the cluster in Channel 2
            cluster_image_channel2 = cluster_mask .* channel2;  
            
            % Combine with Channel 1 and set unused color channels to zero
            cluster_image_combined = cat(3, uint8(channel1_dimmer * 255), uint8(cluster_image_channel2), zeros(size(channel1), 'uint8'));
            
            % Generate k random rotations between 0 and 360 degrees
            for k = 1:6
                % Generate a random angle between 0 and 360 degrees
                random_angle = rand() * 360;
                rotated_cluster_image_combined = imrotate(cluster_image_combined, random_angle, 'crop');
                imshow(rotated_cluster_image_combined);
                
                % Check if any part of the cluster in Channel 2 is out of image range
                rotated_channel2 = rotated_cluster_image_combined(:, :, 2);
                [rotated_height, rotated_width] = size(rotated_channel2);
                cluster_mask_rotated = rotated_channel2 > 0; % Binary mask of the rotated cluster
                if any(cluster_mask_rotated(:))
                    cluster_bounds = regionprops(cluster_mask_rotated, 'BoundingBox');
                    bbox = cluster_bounds(1).BoundingBox;
                    
                    % Check if the cluster bounding box is within image boundaries
                    if bbox(1) < 1 || bbox(2) < 1 || (bbox(1) + bbox(3)) > rotated_width || (bbox(2) + bbox(4)) > rotated_height
                        fprintf('Cluster %d is out of image bounds after rotation %d. Skipping.\n', j, k);
                        continue; % Skip this cluster if it is out of bounds
                    end
                end
                
                % Recalculate the centroid after rotation
                rotated_cluster_center = regionprops(rotated_channel2 > 0, 'Centroid');
                
                if ~isempty(rotated_cluster_center)
                    new_centroid = rotated_cluster_center(1).Centroid;
                    
                    % Determine the translation needed to move the new centroid to the composite center
                    composite_center = [rotated_width / 2, rotated_height / 2]; % Center of the rotated image
                    
                    % Calculate the translation distance
                    translation_x = composite_center(1) - new_centroid(1);
                    translation_y = composite_center(2) - new_centroid(2);
                    
                    % Apply translation using imtranslate with 'OutputView', 'same'
                    aligned_image = imtranslate(rotated_cluster_image_combined, [translation_x, translation_y], 'OutputView', 'same');
                    
                    % Store the aligned rotated cluster image
                    all_rotated_images{end + 1} = aligned_image;
                else
                    fprintf('No valid centroid found for cluster %d after rotation %d.\n', j, k);
                end
            end
        else
            fprintf('Cluster %d does not meet the size criteria. Area: %f\n', j, cluster_properties(j).Area);
        end
    end
end

% Check if any rotated cluster images were processed
if isempty(all_rotated_images)
    error('No rotated cluster images were processed. Please check the input image and the processing criteria.');
end

% Determine the largest size for the final composite image
max_height = 0;
max_width = 0;

% Loop through all images to find the maximum dimensions
for k = 1:length(all_rotated_images)
    [height, width, ~] = size(all_rotated_images{k});
    if height > max_height
        max_height = height;
    end
    if width > max_width
        max_width = width;
    end
end

% Now define the composite image size based on the largest dimensions
composite_image_size = [max_height * 2, max_width * 2, 3]; % Double size for accommodating translations
composite_image = zeros(composite_image_size, 'uint16'); % Use uint16 to avoid overflow when incrementing

% Place each rotated image in the composite image
for k = 1:length(all_rotated_images)
    aligned_image = all_rotated_images{k};
    [height, width, ~] = size(aligned_image);
    
    % Calculate the start position to place the image at the center of the composite
    start_x = max(1, round((composite_image_size(2) - width) / 2));
    start_y = max(1, round((composite_image_size(1) - height) / 2));
    
    % Ensure end indices do not exceed the composite image dimensions
    end_x = min(start_x + width - 1, composite_image_size(2));
    end_y = min(start_y + height - 1, composite_image_size(1));
    
    % Increment the composite image by the original intensities of the aligned images
    composite_image(start_y:end_y, start_x:end_x, 1) = ...
        composite_image(start_y:end_y, start_x:end_x, 1) + uint16(aligned_image(1:(end_y-start_y+1), 1:(end_x-start_x+1), 1));
    composite_image(start_y:end_y, start_x:end_x, 2) = ...
        composite_image(start_y:end_y, start_x:end_x, 2) + uint16(aligned_image(1:(end_y-start_y+1), 1:(end_x-start_x+1), 2));
end

% Normalize the composite image for visualization
composite_image_normalized = uint8(255 * mat2gray(composite_image)) * (1 / 2); % Scale to [0, 255] for display purposes

% Display and save the final merged image
imshow(composite_image_normalized);
imwrite(composite_image_normalized, fullfile(saving_path, 'final_rotate_composite_image.tif'));


