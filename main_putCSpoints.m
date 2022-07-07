%%
close all;
clear ;
clc;

%% load MASK
% Choose result folder
% folderworkspace = pwd;
% folderdestination = uigetdir(folderworkspace,'Specify Results Folder');
% if folderdestination == 0
%   return;
% end

% % Choose preparation folder
% PathName = uigetfile('*.MSK','Specify Prepare Folder');
% if PathName == 0
%     return;
% end
% listing = dir(PathName);
% name = {listing(1:end).name};
% 
% filename_MSK = name(contains(name,'.MSK')==1);
% 
% if length(filename_MSK) == 1
%     target_filename_MSK = fullfile(PathName,filename_MSK);
%     
% else
%     return
% end

% Get directory and file name
PathName = uigetdir(pwd, 'Select Mask Path');
listing = dir(PathName);
name = {listing(1:end).name};
filename_MSK = name(contains(name,'.MSK')==1);
target_filename_MSK = fullfile(PathName,filename_MSK);

% filename_MSK = dir(fullfile(PathName, '*.MSK'));
% filename_MSK = {filename_MSK.name};
mkdir(PathName,'CSData');
folderdestination = fullfile(PathName,'CSData');

for iMSK = 1:size(filename_MSK,2)
%% 計時
t0 = clock;
%% Set intial parameter
main_trans = 0; % main feature transmittance 
main_phas = 0; % main feature phase
main_group = 0; % main feature group(0)
main_poly_info = [main_trans main_phas main_group];
pixel_size_mask = 20; % nano meter
target_pixel_blocks = 5; %5 %10 ;%the number of pixel_size_mask
min_step = 20; % 放點用的最小間距
interval_blocks = 3;
target_pixel_intervals = interval_blocks*target_pixel_blocks; %20;% the number of pixel_size_mask
K_KNN = 5; % K_Block=1


%% cut mask
temp_path = pwd;
% cd(PathName);
% [mask_dimension,sim_region,back_ground_info,msk_poly_target,msk_poly_info,~] = Parse_MSKfile(filename_MSK{iMSK}); % parse MSK information
[mask_dimension,sim_region,back_ground_info,msk_poly_target,msk_poly_info,~] = Parse_MSKfile(target_filename_MSK{1,iMSK}); % parse MSK information

% cd(temp_path);
n_poly_target = length(msk_poly_target);

[Main_Poly,Main_Poly_Info,~] = Main_Polygon(msk_poly_target,msk_poly_info,n_poly_target,sim_region);
[Main_Poly] = mergepolygon(Main_Poly);
[Target_45degree ,Target_90degree ,Target_135degree] = find_45degree(Main_Poly);
n_main_poly_target = length(Main_Poly);
max_poly_decimal = 0;
% [Target_45degree] = find_45degree_merge(Target_45degree,Main_Poly);
[targetCD] = find_CD(Main_Poly);

% Determine the mask grid size
for i = 1:n_main_poly_target
    [row,col] = size(Main_Poly{i});
    for j = 1:row
        for k = 1:col
            poly_str = num2str(Main_Poly{i}(j,k));
            poly_str_dot = find(poly_str == '.');
            if ~isempty(poly_str_dot)
                poly_str_decimal = length(poly_str)-poly_str_dot;
                if  poly_str_decimal > max_poly_decimal
                    max_poly_decimal = poly_str_decimal;
                end
            end
            
        end
    end
end


[target_pixelbased,center_pixel,poly_pixel,meshX,meshY,mask_pixelsizex,mask_pixelsizey,target_pixelbased_blocks,target_pixelbased_intervals,grid_number] = Find_Mask_Center_Pixel2(back_ground_info,sim_region,pixel_size_mask,msk_poly_target,n_poly_target,main_poly_info,target_pixel_blocks,target_pixel_intervals);% polygon to pixel-based

target_pixelnumber = size(target_pixelbased); 

%KNN分類
% [Mask_pixelbased_blocks_classification] = classificationblocks(Mask_pixelbased_blocks);
% [Mask_pixelbased_intervals_classification] = classificationintervals(Mask_pixelbased_intervals);
[Mask_pixelbased_blocks_classification] = classificationblocks(target_pixelbased_blocks,1);
[Mask_pixelbased_intervals_classification] = classificationintervals(target_pixelbased_intervals,K_KNN);

% [Mask_pixelbased_blocks_classification_txt] = classificationblocksusefortxt(target_pixelbased_blocks);

%把intervals的大小對應到blocks的大小
[intervals_1D,intervals_2D] = Resize_Mask_pixelbased_intervals_classification(Mask_pixelbased_intervals_classification,grid_number,interval_blocks);

TargetData.pixel_size_mask = pixel_size_mask;
TargetData.mask_dimension = mask_dimension;
TargetData.sim_region = sim_region;
TargetData.back_ground_info = back_ground_info;
TargetData.msk_poly_target = msk_poly_target;
TargetData.msk_poly_info = msk_poly_info;
TargetData.main_poly_info = main_poly_info;
TargetData.main_poly = Main_Poly;
TargetData.n_poly_target = n_poly_target;
TargetData.n_main_poly_target = n_main_poly_target;
TargetData.target_pixelbased = target_pixelbased;
TargetData.center_pixel = center_pixel; % center of polygon
TargetData.poly_pixel = poly_pixel;
TargetData.meshX = meshX;
TargetData.meshY = meshY;
TargetData.mask_pixelsizex = mask_pixelsizex; % actual pixel size
TargetData.mask_pixelsizey = mask_pixelsizey; % actual pixel size
TargetData.target_pixelnumber = target_pixelnumber; 

% imagesc(abs(TargetData.target_pixelbased));
% B = bwmorph(TargetData.target_pixelbased,'remove',8);
% imshow(B);
% C = TargetData.target_pixelbased(:);
% A = reshape(C,2,2);
% imshow(abs(TargetData.target_pixelbased));
% hold on
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
% end
% L = image()
% mask = boundarymask(L);
% Point_Segment = 25; %10; % in order to have point on polygon
Point_Interval = pixel_size_mask*target_pixel_blocks; %20; % Segment size
% [Polygon_Boundary_Point,Edge_Point] = Find_Polygon_Point(TargetData.n_main_poly_target,TargetData.main_poly,TargetData.sim_region,Point_Segment,Point_Interval); % find the point on the polygon edge
[Polygon_Boundary_Point,polygoncrop_temp,corner_position] = Find_Polygon_Point_v3(TargetData.n_main_poly_target,TargetData.main_poly,TargetData.sim_region,Point_Interval,max_poly_decimal);
%計時
t1 = etime(clock,t0);

[vertexblocks] = Find_vertex_blocks(sim_region,grid_number,Main_Poly);
% [vertexblocks] = Find_vertex_blocks_temp(sim_region,grid_number,Main_Poly);

%CS point資料庫
rule_putCSpoints = [];
[rule_putCSpoints] = ruleputCSpoints(rule_putCSpoints);

%放置CS point
[segmentdata,Polygon_1D_draw,Step] = putCSpoints(Polygon_Boundary_Point,intervals_1D,Mask_pixelbased_blocks_classification,min_step,rule_putCSpoints);

%拿掉polygon裡面的點
[segmentdata] = Takeout_points_Polygon(TargetData.n_main_poly_target,TargetData.main_poly,segmentdata);

% Apply corner radius
CornerRadius = targetCD; %targetCD/2;
% CornerRadius = 125;
Delete_Region = CornerRadius*1.5;
[segmentdata_tmep] = add_corner_radius(TargetData, CornerRadius, segmentdata ,Delete_Region ,Step ,min_step ,vertexblocks ,Target_45degree ,Target_90degree, Target_135degree);

% figure;
% drawPolygon(segmentdata_tmep.cspoints,'*');
% Delete boundary CS point
[segmentdata] = del_bound_pt(segmentdata_tmep,TargetData.sim_region);

%計時
t2 = etime(clock,t0);
t = 'Spend time : ';
disp(t)
disp(t2)


%畫圖
figure;
npoly = Polygon_1D_draw.';
color = [0.75 0.75 0.75];
for i = 1:length(npoly)
    v=npoly{i,1};
    f = 1:length(npoly{i,1});
    patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor',color);
    hold on
end
% drawPolygon(segmentdata.cspoints,'*');
for i = 1: length(segmentdata.cspoints)
scatter(segmentdata.cspoints{i, 1}(:,1),segmentdata.cspoints{i, 1}(:,2),'filled')
hold on
end
% drawPolygon(Polygon_1D_draw);

axis([sim_region(4), sim_region(2), sim_region(3), sim_region(1)]);
title('CS Point')
xlabel('nm')
ylabel('nm')
%% 測試用
% [A] = fragmentation(Polygon_Boundary_Point.poly(17),20,1); %要放有polygon的位置
% drawPolygon(A.cspoints,'*'); %drawPolygon(segmentdata.cspoints,'*');
% drawPolygon(Polygon_Boundary_Point.poly(17));
%% merge MASK
% freeformmask = writemskfile(TargetData,namemask,folderdestination);

%% clipdata established
cssamplelength = 40;
segmentsize = 100;
CSdata = add_csdata(TargetData,segmentsize,cssamplelength,polygoncrop_temp,segmentdata);

[~, current_filename, ~] = fileparts(fullfile(folderdestination, filename_MSK{iMSK}));
save(fullfile(folderdestination,[current_filename '_Step' num2str(min_step) '_K' num2str(K_KNN) '.mat']),'-struct','CSdata');
% clipsData.CSdata = CSdata;
% clipsData.filename_clip = filename_clip;
% [~, filenamesave_clip, ~] = fileparts(filename_clip{1});
% filenamesave_clip = {filenamesave_clip};
% clipsData.filenamesave_clip = filenamesave_clip;
% clipsData.main_poly = TargetData.main_poly;
% clipsData.n_main = length(TargetData.main_poly);
% clipsData.main_center = TargetData.main_center;
% clipsData.n_clip = length(filename_clip);
end