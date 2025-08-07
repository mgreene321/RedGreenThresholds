function [x_target, y_target] = transformImageCoordinates(x_orig, y_orig,dailyToMaster, chooseFiles, dailyCorrPts, masterCorrPts)

if nargin < 6
    transferFolder = uigetdir([], 'Select correspondence points folder');
addpath(genpath(transferFolder));
    if chooseFiles
        if dailyToMaster
            file1 = uigetfile(transferFolder, 'Select daily points');
            file2 = uigetfile(transferFolder, 'Select master points');
        else
            file1 = uigetfile(transferFolder, 'Select master points');
            file2 = uigetfile(transferFolder, 'Select daily points');
            
        end
        
        c1 = importdata(file1);
        c2 = importdata(file2);
        
    else
        if dailyToMaster
            file1 = 'daily_correspondence_points.txt';
            file2 = 'master_correspondence_points.txt';
        else
            file1 = 'master_correspondence_points.txt';
            file2 = 'daily_correspondence_points.txt';
            
        end
        
        fileID = fopen(file1);
        c1 = fscanf(fileID, '%f, %f', [2 Inf]);
        c1 = c1';
        fclose(fileID);
        
        fileID = fopen(file2);
        c2 = fscanf(fileID, '%f, %f', [2 Inf]);
        c2 = c2';
        fclose(fileID);
    end
else % input daily and master correspondence points
    
    if dailyToMaster
    c1 = dailyCorrPts;
    c2 = masterCorrPts;
    else
        c2 = dailyCorrPts;
    c1 = masterCorrPts;
        
    end
    
    
end

for i = 1:numel(x_orig)
    try
        [x_target(i,:), y_target(i,:)] = transformSingleCoordinate(x_orig(i,:), y_orig(i,:), c1, c2);
    catch
        x_target(i,:) = nan; y_target(i,:) = nan;
    end
end



function [x2, y2] = transformSingleCoordinate(x1, y1,c1,c2)

DT1 = delaunayTriangulation(c1);
%DT2 = delaunayTriangulation(c2);

triIdx = DT1.pointLocation(x1, y1); %get coordinates of vertices that enclose point

rowIndices = DT1(triIdx,:);

A = computeAffine(c1(rowIndices,:), c2(rowIndices,:));
temp = A*[x1 y1 1]';
x2 = temp(1);
y2 = temp(2);

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
