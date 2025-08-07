function [Raw_image_fundus_scaled, coneData] = prepareOCTImageForConeTransfer(subjectId, dailyPPD)

%tempDir = dir(fullfile(pathToConeIDMat, 'Cone_identification*.mat'));
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
%subjectFolder = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\MissingConeIdentificationData\20210R';

OCTimDir = dir(fullfile(subjectFolder, 'OCTMaster', '*Cone_identification*.mat'));

% Put OCT image into fundus view orientation
load(fullfile(OCTimDir(1).folder, OCTimDir(1).name));
Raw_image_fundus= fliplr(rot90(Raw_image,-1));
masterPPD = round(mean(size(Raw_image)./FOV_in_degrees));
scaleFactor = dailyPPD/masterPPD;

% Scale OCT image 

numRows1 = size(Raw_image_fundus,1);
numCols1 = size(Raw_image_fundus,2);

numRows2 = ceil(scaleFactor*numRows1);
numCols2 = ceil(scaleFactor*numCols1);

[X1, Y1] = meshgrid(1:numCols1, 1:numRows1);
[X2, Y2] = meshgrid(1:numCols2, 1:numRows2);

Xgrid= X1 .* scaleFactor;
Ygrid = Y1 .* scaleFactor;

Raw_image_fundus_scaled = interp2(Xgrid, Ygrid, Raw_image_fundus, X2, Y2, 'cubic');


%transform cone coordinates

%Rotation matrix
theta = pi/2;
R = [cos(theta) -sin(theta)
    sin(theta) cos(theta)];

cone_locs = str2double(Cone_Mat_all_marked(:, [1 2]));

% Transform every location in cone_locs with the transformation matrix
for i = 1:size(Cone_Mat_all_marked,1)
    x = cone_locs(i,1);
    y = cone_locs(i,2);
    xy1 = [x y];
    xy2(i,:) = transpose(R*xy1');
    
end

% translate

xy2(:,1) = xy2(:,1) - (min(xy2(:,1))-min(cone_locs(:,2)));
xy2(:,2) = xy2(:,2) - (min(xy2(:,2))-min(cone_locs(:,1)));

% reflect across vertical dimension

xy2(:,1) = -xy2(:,1);

%translate 
xy2(:,1) = xy2(:,1) - (min(xy2(:,1))-min(cone_locs(:,2)));

% scale

x2grid = scaleFactor.*(1:numCols1);
y2grid = scaleFactor.*(1:numRows1);

x2 = interp1(1:numCols1, x2grid, xy2(:,1), 'cubic');
y2 = interp1(1:numRows1, y2grid, xy2(:,2), 'cubic');


coneData = [x2 y2 Cone_Mat_all_marked(:,3)];

savePath = fullfile(subjectFolder, 'OCTMaster');
%savePath = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\MissingConeIdentificationData\20210R\OCTMaster';
imwrite(Raw_image_fundus_scaled./max(Raw_image_fundus_scaled(:)), fullfile(savePath, 'Raw_image_fundus_scaled.tif'));

fileName = fullfile(savePath, 'coneData.mat');


save(fileName, 'coneData', 'Raw_image_fundus_scaled');


end
