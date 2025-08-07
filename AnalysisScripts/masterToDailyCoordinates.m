function [x_daily, y_daily]=masterToDailyCoordinates(x_master,y_master, subjectId, expDate, plotFlag)
root = getenv('USERPROFILE');

%tempDir = dir(fullfile(root, 'Documents', '**', '*RedGreenThresholds', 'ConeTransfer', subjectId));

tempDir = dir(fullfile(root, 'Documents', '**', '*RedGreenThresholds', subjectId, expDate, 'ConeTransfer'));
transferFullPath = tempDir(1).folder;

% if ~exist('transferDataFolder', 'var') || isempty(transferDataFolder)
%     transferDataFolder = uigetdir(coneTransferRoot);
% end


if ~exist('plotFlag', 'var') || isempty(plotFlag)
    plotFlag = 0;
end
% On PC

if plotFlag
    h = figure; hold on;
    masterDir = dir([transferFullPath filesep '*croppedMaster*']);
    dailyDir = dir([transferFullPath filesep '*customRef*']);
    
    masterIm = imread(fullfile(masterDir.folder,masterDir.name));
    dailyIm = imread(fullfile(dailyDir.folder,dailyDir.name));
    
    
    set(h, 'Position', [h.Position(1) h.Position(2) 2.*h.Position(3) h.Position(4)]);
    
    subplot(1, 2, 2); imagesc(dailyIm), colormap gray; hold on; axis off;
    subplot(1, 2, 1); imagesc(masterIm), colormap gray; hold on; axis off;
end
dailyFullFile = fullfile(transferFullPath, 'daily_correspondence_points.txt');
masterFullFile = fullfile(transferFullPath, 'master_correspondence_points.txt');
rectDir = dir(fullfile(transferFullPath, 'rect*'));
rect = importdata(fullfile(transferFullPath,rectDir.name));
rect = [0 0 0 0];

[dc, mc] = getCorrespondencePts(dailyFullFile, masterFullFile);
DT_daily = delaunayTriangulation(dc);
DT_master = delaunayTriangulation(mc);


x_daily = nan(size(x_master)); y_daily = nan(size(y_master));
for i = 1:length(x_master)
    
    try
    %x_master = x_master - rect(1)+1; y_master = y_master - rect(2)+1;
    x = x_master(i) - rect(1)+1; y = y_master(i) - rect(2)+1;
    
    triIdx = DT_master.pointLocation(x, y); %get coordinates of vertices that enclose point
    rowIndices = DT_master(triIdx,:);
    
    A = computeAffine(mc(rowIndices,:), dc(rowIndices,:));
    temp = A*[x y 1]';
    x_daily(i) = temp(1); y_daily(i) = temp(2);
    %x_daily = temp(1) + rect(1)-1; y_daily = temp(2) + rect(2)-1;
    
    catch
    end
    
    if plotFlag
        
        
        subplot(1, 2, 2),scatter(x_daily(i), y_daily(i), 'yo'); title('Daily Map');
        
        
        subplot(1,2,1), scatter(x, y, 'yo'); title('Master Map');
        
    else
    end
    
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