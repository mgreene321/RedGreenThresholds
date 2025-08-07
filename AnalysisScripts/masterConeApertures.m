function [coneApertures, cone_locs, cone_labels] = masterConeApertures(subjectId)
coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 11;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma);
root = getenv('USERPROFILE');

%tempDir = dir(fullfile(root, 'Documents', '**', '*RedGreenThresholds', subjectId, expDate, 'ConeTransfer'));
tempDir = dir(fullfile(root, 'Documents', '**', '*RedGreenThresholds', subjectId, 'OCTMaster'));

%coneTransferRoot = tempDir(1).folder;
%selpath = uigetdir(coneTransferRoot); % Select folder with relevant cone transfer data
transferFullPath = tempDir(1).folder;
%dailyDir = dir([transferFullPath filesep '*customRef*']);
retina_map = imread(fullfile(transferFullPath,'Raw_image_fundus_scaled.tif'));

load(fullfile(transferFullPath,'coneData.mat'));
cone_locs = str2double(coneData(:, [1 2]));
cone_labels = coneData(:,end);

nL = sum(cone_labels(~any(isnan(cone_locs),2)) == "L");
nM = sum(cone_labels(~any(isnan(cone_locs),2)) == "M");
nS = sum(cone_labels(~any(isnan(cone_locs),2)) == "S");

N = nL + nM + nS;
pL = nL./N;
pM = nM./N;
pS = nS./N;

load(fullfile(transferFullPath, 'likelyMissingCones.mat'));

% Get master cone locs

%
% jsonDir = dir([transferFullPath '\*2k2k_region.json']);
% jsonfile = fullfile(jsonDir.folder,jsonDir.name);
% str = fileread(jsonfile);
% data = jsondecode(str);
% cone_data = transpose([data.cone_data.l]);
% cone_locs = cone_data(:, [1 2]); %N x 2 matrix, rows are cones, columns x and y coords
% cone_labels = [data.cone_data.c]'; %N x 1 matrix with cone types

L_apertures = zeros(size(retina_map));
M_apertures = zeros(size(retina_map));
S_apertures = zeros(size(retina_map));

for i = 1:size(cone_locs,1) % for each cone
    Xc = round(cone_locs(i,1)); Yc = round(cone_locs(i,2));
    
    if Xc - fix(coneSpacingPixels/2) > 0 && Yc - fix(coneSpacingPixels/2) > 0
        
        if Xc + fix(coneSpacingPixels/2) < size(retina_map,2) &&...
                Yc + fix(coneSpacingPixels/2) < size(retina_map,1)
            
            if strcmpi(cone_labels(i), 'l') %L submosaic
                L_temp  = zeros(size(retina_map));
                L_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = coneAperture;
                L_apertures = max(cat(3, L_apertures, L_temp),[],3);
                %L_apertures = L_apertures + L_temp;
            elseif strcmpi(cone_labels(i), 'm') %M submosaic
                M_temp  = zeros(size(retina_map));
                M_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = coneAperture;
                M_apertures = max(cat(3, M_apertures, M_temp),[],3);
                %M_apertures = M_apertures + M_temp;
                
            elseif strcmpi(cone_labels(i), 's') %S submosaic
                S_temp  = zeros(size(retina_map));
                S_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = coneAperture;
                S_apertures = max(cat(3, S_apertures, S_temp),[],3);
                %S_apertures = S_apertures + S_temp;
            else %missing
                L_temp  = zeros(size(retina_map));
                L_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pL.*coneAperture;
                L_apertures = max(cat(3, L_apertures, L_temp),[],3);
                M_temp  = zeros(size(retina_map));
                M_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pM.*coneAperture;
                M_apertures = max(cat(3, M_apertures, M_temp),[],3);
                S_temp  = zeros(size(retina_map));
                S_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pS.*coneAperture;
                S_apertures = max(cat(3, S_apertures, S_temp),[],3);
                
                
            end
        else
        end
    end
end

for i = 1:size(likelyMissingCones,1)
     Xc = round(likelyMissingCones(i,1)); Yc = round(likelyMissingCones(i,2));
    
     if Xc - fix(coneSpacingPixels/2) > 0 && Yc - fix(coneSpacingPixels/2) > 0
         try
    L_temp  = zeros(size(retina_map));
    L_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pL.*coneAperture;
    L_apertures = max(cat(3, L_apertures, L_temp),[],3);
    M_temp  = zeros(size(retina_map));
    M_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pM.*coneAperture;
    M_apertures = max(cat(3, M_apertures, M_temp),[],3);
    S_temp  = zeros(size(retina_map));
    S_temp(Yc-fix(windowSize/2):Yc+fix(windowSize/2), Xc-fix(windowSize/2):Xc+fix(windowSize/2)) = pS.*coneAperture;
    S_apertures = max(cat(3, S_apertures, S_temp),[],3);
         catch
         end
     end
end

coneApertures(:,:,1) = L_apertures;
coneApertures(:,:,2) = M_apertures;
coneApertures(:,:,3) = S_apertures;

coneApertures(coneApertures > max(coneAperture(:))) = max(coneAperture(:));
save(fullfile(transferFullPath, 'coneApertures.mat'), 'coneApertures');
end