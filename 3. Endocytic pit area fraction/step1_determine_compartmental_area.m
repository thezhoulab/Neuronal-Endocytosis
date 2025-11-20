%% This file is used to identify cell masks for compartmental area

clear all
close all

parent_folder = 'Y:\20241208 D54D2 in b2KD and WT_APP\Analysis';
dirs = dir(parent_folder);

for m = 1%:length(dirs)
    if dirs(m).isdir && ~strcmp(dirs(m).name, '.') && ~strcmp(dirs(m).name, '..')
        folder_path = fullfile(parent_folder, dirs(m).name);
        saving_path = fullfile(folder_path, 'step 1_cell mask_MAP2');
        mkdir(saving_path);

        file_list = dir(fullfile(folder_path, 'MAP2_AF647*.tif')); % adjust the name

        for i = 1:length(file_list)
            img_raw = imread(fullfile(folder_path, file_list(i).name));
            img = im2double(img_raw);

            filter = fspecial('average', 3);
            img = imfilter(img, filter, 'replicate');
            img = imbinarize(img, "adaptive", 'sensitivity', 0.4); % adjust sensitivity

            boundaries = bwboundaries(img);
            if isempty(boundaries)
                warning('No boundary detected in %s', file_list(i).name);
                continue
            end

            % Find largest boundary
            largestBoundaryIndex = 0;
            largestBoundaryArea = 0;
            for k = 1:numel(boundaries)
                boundary = boundaries{k};
                boundaryArea = size(boundary, 1);
                if boundaryArea > largestBoundaryArea
                    largestBoundaryIndex = k;
                    largestBoundaryArea = boundaryArea;
                end
            end

            largestBoundary = boundaries{largestBoundaryIndex};
            mask = poly2mask(largestBoundary(:,2), largestBoundary(:,1), size(img,1), size(img,2));
            mask = imfill(mask, 'holes');
            img = img .* mask;

            se = strel('disk', 5);
            img = imclose(img, se);
            cc = bwconncomp(img);
            min_size = 10000; % adjust min_size
            for obj = 1:cc.NumObjects
                if numel(cc.PixelIdxList{obj}) < min_size
                    img(cc.PixelIdxList{obj}) = 0;
                end
            end

            figure
            subplot(1,2,1), imshow(imadjust(img_raw)), title('original');
            subplot(1,2,2), imshow(img), title('cell mask');

            filename = fullfile(saving_path, sprintf('cellmask_%04d.tif', i));
            imwrite(img, filename, 'Compression', 'none');
        end
    end
end
