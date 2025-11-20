%% This file is used to differentiate surface-localized and internalized endocytic pits
%Type 1 is surface-localized
%Type 2 is internalized

clear all;
close all;
clc;

% Parameters for rendering and DBSCAN
epsilon = 0.25; % Maximum distance between points (nm)
minPts = 10; % Minimum number of points per cluster
width = 0.1; % Gaussian width for rendering
a1 = 20; % Image scale for Channel 1
a2 = 20; % Image scale for Channel 2
normalization_factor = 0.5; % Normalize z-coordinates

folder_path = '\\E1-054867\Jinyu\20250927 DIV14 OE_EndoA2_583 b2_647\Analysis\diat_regions'; 
filename = dir(fullfile(strcat([folder_path '\'], '*.bin')));
saving_path = fullfile(folder_path, '2 Type Puncta');
mkdir(saving_path);

% Initialize result storage
results_summary = fopen(fullfile(saving_path, 'results_summary.txt'), 'w');
fprintf(results_summary, 'File Name\tSeg_Num(#)\tType 1(#)\tType 2(#)\n');

for ii = 1:length(filename)
    % Import and read molecule list
    try
        [MList, memoryMap] = ReadMasterMoleculeList([folder_path '\' filename(ii).name], ...
            'fieldsToLoad', {'xc', 'yc', 'zc', 'c'}, 'ZScale', 160);
    catch
        fprintf('Error reading file: %s. Skipping...\n', filename(ii).name);
        continue;
    end
    
    % Extract coordinates for each channel
    xcdata_1 = MList.xc(MList.c == 1);
    ycdata_1 = MList.yc(MList.c == 1);
    zcdata_1 = MList.zc(MList.c == 1);
    
    xcdata_2 = MList.xc(MList.c == 2);
    ycdata_2 = MList.yc(MList.c == 2);
    zcdata_2 = MList.zc(MList.c == 2);
    
    % Skip if no data exists
    if isempty(xcdata_1) || isempty(xcdata_2)
        fprintf('No data in file %s, skipping...\n', filename(ii).name);
        continue;
    end
    
    % Segment the dataset into 3 pieces based on xcdata_1
    xcdata1_min = min(xcdata_1);
    xcdata1_max = max(xcdata_1);
    num_segments = 1; % 1 for all
    xcdata1_edges = linspace(xcdata1_min, xcdata1_max, num_segments + 1);
    
    for seg = 1:num_segments
        % Define segment boundaries
        seg_min = xcdata1_edges(seg);
        seg_max = xcdata1_edges(seg + 1);
        
        % Filter data for Channel 1
        seg_idx_1 = xcdata_1 >= seg_min & xcdata_1 < seg_max;
        xcdata_1_seg = xcdata_1(seg_idx_1);
        ycdata_1_seg = ycdata_1(seg_idx_1);
        zcdata_1_seg = zcdata_1(seg_idx_1);
        
        % Filter data for Channel 2
        seg_idx_2 = xcdata_2 >= seg_min & xcdata_2 < seg_max;
        xcdata_2_seg = xcdata_2(seg_idx_2);
        ycdata_2_seg = ycdata_2(seg_idx_2);
        zcdata_2_seg = zcdata_2(seg_idx_2);
        
        % Skip if no data exists in this segment
        if length(xcdata_1_seg) < 10 || length(xcdata_2_seg) < 10
            fprintf('Not enough data in segment %d, skipping...\n', seg);
            continue;
        end
        
        % Normalize Z-coordinates
        zcdata_1_normalized = zcdata_1_seg / 1200 * range(xcdata_1) * normalization_factor;
        zcdata_2_normalized = zcdata_2_seg / 1200 * range(xcdata_1) * normalization_factor;
        
        % 3D Puncta Identification
        data_3D = [xcdata_2_seg, ycdata_2_seg, zcdata_2_normalized];
        labels = dbscan(data_3D, epsilon, minPts);
        
        % Initialize stacks and counters
        puncta_membrane_stack = [];
        moleculeCountType1 = 0;
        moleculeCountType2 = 0;

        % Process each puncta
        unique_clusters = unique(labels(labels > 0)); % Exclude noise
        for cluster_id = unique_clusters'
            % Get the current puncta
            puncta_points_idx = find(labels == cluster_id);
            puncta_xc = xcdata_2_seg(puncta_points_idx);
            puncta_yc = ycdata_2_seg(puncta_points_idx);
            puncta_zc = zcdata_2_seg(puncta_points_idx);
            puncta_zc = puncta_zc / 1200 * range(xcdata_1) * normalization_factor;
            
            % Calculate the convex hull volume for the current puncta
            puncta_coords = [puncta_xc, puncta_yc, puncta_zc]; % 3D coordinates of the puncta
            if size(puncta_coords, 1) >= 4 % Convex hull requires at least 4 points
                [~, volume] = convhull(double(puncta_coords));
            else
                volume = 0; % If fewer than 4 points, volume cannot be calculated
            end
            
            % Create a YZ projection for the current puncta
            ROI_yz = [min(zcdata_1_normalized), max(zcdata_1_normalized); ...
                min(ycdata_1_seg), max(ycdata_1_seg)];
            renderedStack_yz_1 = RenderMList([ycdata_1_seg, zcdata_1_normalized], ...
                'gaussianWidth', width, 'ROI', ROI_yz, 'imageScale', a1);
            renderedStack_yz_puncta = RenderMList([puncta_yc, puncta_zc], ...
                'gaussianWidth', width, 'ROI', ROI_yz, 'imageScale', a2);
            
            % Calculate overlay area
            bw_channel_1 = imbinarize(renderedStack_yz_1);
            
            % Ellipse Fitting
            bw_channel_1 = imfill(bw_channel_1, 'holes'); % Fill any holes
            % Extract boundaries of all connected components
            all_boundaries = bwboundaries(bw_channel_1);
            
            % Combine boundary points from all components
            combined_boundary_points = [];
            for k = 1:length(all_boundaries)
                combined_boundary_points = [combined_boundary_points; all_boundaries{k}];
            end
            
            % Check if there are enough boundary points
            if size(combined_boundary_points, 1) < 5
                fprintf('Too few boundary points for segment %d in file %s. Skipping...\n', seg, filename(ii).name);
                continue;
            end
            
            % Extract Y and Z coordinates of combined boundary points
            y_boundary = combined_boundary_points(:, 1); % Y-coordinates
            z_boundary = combined_boundary_points(:, 2); % Z-coordinates
            
            % Fit an ellipse to the combined boundary points using `fit_ellipse`
            ellipse_t = fit_ellipse(z_boundary, y_boundary);
            if isempty(ellipse_t)
                fprintf('Ellipse fitting failed for segment %d in file %s. Skipping...\n', seg, filename(ii).name);
                continue;
            end
            
            % Identify Puncta in Channel 2
            bw_channel_2 = imbinarize(renderedStack_yz_puncta);
            
            % Measure properties of connected components in the binary mask
            props = regionprops(bw_channel_2, 'Area', 'Centroid', 'PixelIdxList');
            
            % Check if there are any components
            if isempty(props)
                warning('No connected components detected in bw_channel_2!');
                continue; % Return empty if no components found
            else
                % Find the largest connected component by Area
                [~, largest_idx] = max([props.Area]);
                
                % Keep only the largest component
                largest_segment_mask = false(size(bw_channel_2));
                largest_segment_mask(props(largest_idx).PixelIdxList) = true;
                bw_channel_2 = largest_segment_mask;
                
                % Get the centroid of the largest component
                puncta_centroid = props(largest_idx).Centroid; % [Z, Y] coordinates
            end
            
            % Calculate d1 and d2
            [d1, d2] = calculateDistances(puncta_centroid, ellipse_t);
            
            % Classify puncta based on d1 and d2
            if d1 >= (3/8) * d2 && d1 <= (1/2) * d2
                puncta_membrane_stack = [puncta_membrane_stack; cluster_id];
                
                ROI_xy = [min(ycdata_1_seg), max(ycdata_1_seg); ...
                    min(xcdata_1_seg), max(xcdata_1_seg)];
                renderedStack_xy_1 = RenderMList([xcdata_1_seg, ycdata_1_seg], ...
                    'gaussianWidth', width, 'ROI', ROI_xy, 'imageScale', a1);
                renderedStack_xy_puncta = RenderMList([puncta_xc, puncta_yc], ...
                    'gaussianWidth', width, 'ROI', ROI_xy, 'imageScale', a2);
                
                % Binarize XY projections
                imbw_1 = imbinarize(renderedStack_xy_1);
                imbw_2 = imbinarize(renderedStack_xy_puncta);
                imbw_1 = double(imbw_1);
                imbw_2 = double(imbw_2);
                
                % Ensure images are cropped to the same area
                min_length = min(size(imbw_1, 1), size(imbw_2, 1));
                min_width = min(size(imbw_1, 2), size(imbw_2, 2));
                imbw_1 = imbw_1(1:min_length, 1:min_width);
                imbw_2 = imbw_2(1:min_length, 1:min_width);
                
                % Crop to relevant region based on ch1 presence (optional)
                validRows = any(imbw_1, 2) | any(imbw_2, 2);
                imbw_1 = imbw_1(validRows, :);
                imbw_2 = imbw_2(validRows, :);
                
                % Get puncta pixels (assuming a single puncta in imbw_2)
                puncta_pixels = find(imbw_2);
                puncta_area = numel(puncta_pixels); % Number of pixels
                overlap_area = sum(imbw_1(puncta_pixels));
                
                % Classification logic for this single puncta
                moleculeCountType1 = moleculeCountType1 + 1;
                      
            elseif d1 <= (3/8)* d2
                moleculeCountType2 = moleculeCountType2 + 1;
            end
        end
        
    end
    
    % Write to results summary (with file name and segment)
    fprintf(results_summary, '%s\t%d\t%d\t%d\n', ...
        filename(ii).name, seg, moleculeCountType1, moleculeCountType2);
    
    % Visualization: Render Final Planes
    ROI_xy = [min(ycdata_1_seg), max(ycdata_1_seg); min(xcdata_1_seg), max(xcdata_1_seg)];
    ROI_xz = [min(zcdata_1_normalized), max(zcdata_1_normalized); min(xcdata_1_seg), max(xcdata_1_seg)];
    ROI_yz = [min(zcdata_1_normalized), max(zcdata_1_normalized); min(ycdata_1_seg), max(ycdata_1_seg)];
    
    % Render planes for Channel 1
    renderedStack_xy_1 = RenderMList([xcdata_1_seg, ycdata_1_seg], 'gaussianWidth', width, 'ROI', ROI_xy, 'imageScale', a1);

    % Render planes for 'membrane localized' and 'inside localized'
    puncta_xc_localized = xcdata_2_seg(ismember(labels, [puncta_membrane_stack]));
    puncta_yc_localized = ycdata_2_seg(ismember(labels, [puncta_membrane_stack]));
    
    % Render planes for noise puncta
    renderedStack_xy_2_localized = RenderMList([puncta_xc_localized, puncta_yc_localized], 'gaussianWidth', width, 'ROI', ROI_xy, 'imageScale', a2);
    
    % Combine and save planes
    combined_rgb_xy = cat(3, renderedStack_xy_1, renderedStack_xy_2_localized, zeros(size(renderedStack_xy_1)));
    imwrite(uint8(255 * mat2gray(combined_rgb_xy)), fullfile(saving_path, sprintf('%s_segment%d_combined_rgb_xy.tif', filename(ii).name, seg)));
    
end
fclose(results_summary);
disp('Analysis complete.');


% Function: Calculate d1 and d2
function [d1, d2] = calculateDistances(puncta_centroid, ellipse_t)
% Extract ellipse parameters
xc = ellipse_t.X0_in;
yc = ellipse_t.Y0_in;
a = ellipse_t.a;
b = ellipse_t.b;
phi = ellipse_t.phi;

% Transform puncta centroid to rotated ellipse frame
R = [cos(phi), sin(phi); -sin(phi), cos(phi)];
puncta_transformed = R * [(puncta_centroid(2) - xc); (puncta_centroid(1) - yc)];

zp_transformed = puncta_transformed(1);
yp_transformed = puncta_transformed(2);

% Calculate d1
d1 = sqrt(zp_transformed^2 + yp_transformed^2);

% Ellipse intersection points
slope = yp_transformed / zp_transformed;
A = b^2 + a^2 * slope^2;
C = -a^2 * b^2;
z_cross = sqrt(-C / A) * [-1, 1];
y_cross = slope * z_cross;

% Transform intersection points back to original frame
cross_points_transformed = R' * [z_cross; y_cross];
z_cross_original = cross_points_transformed(1, :) + xc;
y_cross_original = cross_points_transformed(2, :) + yc;

% Calculate d2
d2 = sqrt((z_cross_original(1) - z_cross_original(2))^2 + ...
    (y_cross_original(1) - y_cross_original(2))^2);
end
