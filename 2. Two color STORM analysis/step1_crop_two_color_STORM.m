%% This file is used to crop regions for two-color STORM analysis

clear all;
close all;
clc;
binpath = '\\E1-054867\Jinyu\20250324 b3_AF647 CHC_CF583R\Analysis'; %%%%%%% path that you saved your png and bin files
mkdir([binpath '\diat_regions']);
ipath=[binpath '\diat_regions\'];
filename = dir(fullfile(strcat([binpath '\'],'*_crop.png')));%%%% the png must ends with "_crop", otherwise change this line
filename_bin = dir(fullfile(strcat([binpath '\'],'*_list_corr.bin')));

for ii=1:length(filename_bin)
B = regexp(filename_bin(ii).name,'\d*','Match');
for i= 1:length(B)
  if ~isempty(B{i})
      Num(i,1)=str2double(B{i});
  else
      Num(i,1)=NaN;
  end
end

image_file=[num2str(Num(end)) '_crop.png'];    
binname = filename_bin(ii).name;
% if ~exist([binpath image_file], 'file')
%     continue;
% end
% data = importdata([binpath '\' binname]);
[mList, memoryMap] = ReadMasterMoleculeList([binpath '\' binname]);
Mlist=[];
for i=1:length(mList.x)
    %,'fieldsToLoad',{'xc','yc','z'},'ZScale',167);
    Mlist(i).x=mList.x(i);
    Mlist(i).y=mList.y(i);
    Mlist(i).xc=mList.xc(i);
    Mlist(i).yc=mList.yc(i);
    Mlist(i).h=mList.h(i);
    Mlist(i).a=mList.a(i);
    Mlist(i).w=mList.w(i);
    Mlist(i).phi=mList.phi(i);
    Mlist(i).ax=mList.ax(i);
    Mlist(i).bg=mList.bg(i);
    Mlist(i).i=mList.i(i);
    Mlist(i).c=mList.c(i);
    Mlist(i).density=mList.density(i);
    Mlist(i).frame=mList.frame(i);
    Mlist(i).length=mList.length(i);
    Mlist(i).link=mList.link(i);
    Mlist(i).z=mList.z(i);
    Mlist(i).zc=mList.zc(i);
end
% data = importdata([binpath binname(1:end-8) 'drift.txt']);
%   x_drift = data(:,2); %pixel -> nm
%   y_drift = data(:,3);
%   z_drift = data(:,4);
% [x_drift,y_drift] = XcorrDriftCorrect(mList);
% for i=1:length(Mlist)
% Mlist(i).xc = Mlist(i).x - x_drift(Mlist(i).frame);
% Mlist(i).yc = Mlist(i).y - y_drift(Mlist(i).frame);
% Mlist(i).zc = Mlist(i).z - z_drift(Mlist(i).frame);
% end

im=double(imread([binpath '\' image_file]));
[im_width,im_length,color]=size(im);
figure; imshow(im);hold on
crop_image=zeros(im_width,im_length);
[crop_index_row,crop_index_col]=find(im(:,:,1)==0&im(:,:,2)==0&im(:,:,3)==255);% blue
for i=1:length(crop_index_row)
crop_image(crop_index_row(i),crop_index_col(i))=1;
end
crop_image = imfill(crop_image,'holes');
crop_image=bwareaopen(crop_image,100);

figure;
imshow(crop_image);
imbinarize(double(crop_image));

crop_image_label=bwlabel(crop_image);
imshow(crop_image_label);
data_crop=[];
%% 
for j=1:length(Mlist)
Mlist(j).boxnumber=0;
end


for i=1:max(max(crop_image_label))
    [crop_index_row,crop_index_col]=find(crop_image_label==i);
    index=find([Mlist.xc]*im_length/256>min(crop_index_col)&[Mlist.xc]*im_length/256<max(crop_index_col)&[Mlist.yc]*im_length/256<max(crop_index_row)&[Mlist.yc]*im_length/256>min(crop_index_row));
    for j=index
        if ismember([round(Mlist(j).xc*im_length/256),round(Mlist(j).yc*im_length/256)],[crop_index_col,crop_index_row],'rows');
          Mlist(j).boxnumber=i;
         else Mlist(j).boxnumber=0;
        end
    end
end
%% 

for i=1:max(max(crop_image_label))
    data_crop=Mlist(find([Mlist.boxnumber]==i));
     figure; scatter([data_crop.x]*im_length/256,[data_crop.y]*im_length/256,'.');
    [crop_index_row,crop_index_col]=find(crop_image_label==i);
    crop_center_x=mean(crop_index_row);
    crop_center_y=mean(crop_index_col);
     distance=pdist2([crop_center_x,crop_center_y],[crop_index_row,crop_index_col]);
     [sorteddis,sortingindices]=sort(distance, 'descend');
     corner_x=crop_index_row(sortingindices(1:10));
     corner_y=crop_index_col(sortingindices(1:10));

    corner_x=max(crop_index_row);
    corner_y=min(crop_index_col(find(crop_index_row==corner_x)));    
    theta_1=atan((corner_y-crop_center_y)./(corner_x-crop_center_x))/pi*180;
    if (corner_y-crop_center_y)<0 & (corner_x-crop_center_x)<0
        theta_1=theta_1+180;
    end
%      angle=unique(round(theta/20));
%      if length(angle)<2
%          continue
%     end
%      theta_1=mean(theta(find(round(theta/20)==angle(1))));
%      theta_2=mean(theta(find(round(theta/20)==angle(2))));
    corner_y=max(crop_index_col);
    corner_x=max(crop_index_row(find(crop_index_col==corner_y))); 
    theta_2=atan((corner_y-crop_center_y)./(corner_x-crop_center_x))/pi*180;
    
    if abs(theta_1-theta_2)>=90
        rotate_angle=-(theta_1+theta_2)/2;
    else rotate_angle=90-(theta_1+theta_2)/2;
    end
    rotate_angle=-deg2rad(rotate_angle);
    R=[cos(rotate_angle),-sin(rotate_angle);sin(rotate_angle), cos(rotate_angle)];
    coordinate_new=R*[[data_crop.xc];[data_crop.yc]];
    for jj=1:length(coordinate_new)
    data_crop(jj).xc=coordinate_new(1,jj)-min(coordinate_new(1,:))+5;
    data_crop(jj).yc=coordinate_new(2,jj)-min(coordinate_new(2,:))+5;
    end
    WriteMoleculeList(data_crop, [ipath binname '_' num2str(i) '.bin'])
%     clear data_crop
    
end
end