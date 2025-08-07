function [x_master, y_master]=dailyToMasterCoordinates(x_daily,y_daily, plotFlag)

if nargin < 3
    plotFlag = 0;
end

% On PC

root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
selpath = uigetdir(RedGreenThresholdsPath); % Select folder with relevant cone transfer data

dailyFullFile = fullfile(selpath, 'daily_correspondence_points.txt');
masterFullFile = fullfile(selpath, 'master_correspondence_points.txt');

[dc, mc] = getCorrespondencePts(dailyFullFile, masterFullFile);
DT_daily = delaunayTriangulation(dc);
DT_master = delaunayTriangulation(mc);

triIdx = DT_daily.pointLocation(x_daily, y_daily); %get coordinates of vertices that enclose point

rowIndices = DT_daily(triIdx,:);

A = computeAffine(dc(rowIndices,:), mc(rowIndices,:));
temp = A*[x_daily y_daily 1]';
x_master = temp(1); y_master = temp(2);

if plotFlag
    
    %masterDir = dir([selpath filesep '*warped_master_texture*']);
    %dailyDir = dir([selpath filesep 'output_texture*.png']);
    [masterImFile, masterImPath] = uigetfile(selpath, 'Select Master Image');
    [dailyImFile, dailyImPath] = uigetfile(selpath, 'Select Daily Image');
    
    
    masterIm = imread(fullfile(masterImPath, masterImFile));
    dailyIm = imread(fullfile(dailyImPath, dailyImFile));
    
    h = figure; hold on;
    set(h, 'Position', [h.Position(1) h.Position(2) 2.*h.Position(3) h.Position(4)]);
    
    subplot(1, 2, 1); imagesc(dailyIm), colormap gray; hold on; axis off;
    scatter(x_daily, y_daily, 'yo'); title('Daily Map');
    
    subplot(1, 2, 2); imagesc(masterIm), colormap gray; hold on; axis off;
    scatter(temp(1), temp(2), 'yo'); title('Master Map');
    
else
end


function [daily_correspondence, master_correspondence] = ...
    getCorrespondencePts(dailyFullFile, masterFullFile)

fileID = fopen(dailyFullFile);
daily_correspondence = fscanf(fileID, '%f, %f', [2 Inf]);
daily_correspondence = daily_correspondence';
fclose(fileID);

fileID = fopen(masterFullFile);
master_correspondence = fscanf(fileID, '%f, %f', [2 Inf]);
master_correspondence = master_correspondence';
fclose(fileID);

function A = computeAffine(tri1_pts, tri2_pts)

%Written by James Fong <3
origin_to_tri1 = eye(3);
origin_to_tri1(1:2,1) = tri1_pts(2,:) - tri1_pts(1,:);
origin_to_tri1(1:2,2) = tri1_pts(3,:) - tri1_pts(1,:);
origin_to_tri1(1:2,3) = tri1_pts(1,:);

origin_to_tri2 = eye(3);
origin_to_tri2(1:2,1) = tri2_pts(2,:) - tri2_pts(1,:);
origin_to_tri2(1:2,2) = tri2_pts(3,:) - tri2_pts(1,:);
origin_to_tri2(1:2,3) = tri2_pts(1,:);

tri1_to_origin = inv(origin_to_tri1);
A = origin_to_tri2 * tri1_to_origin;