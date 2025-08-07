function [customRef, croppedMasterMap] = prepareImagesForConeTransfer

% The custom reference should be 512 x 512
imFileTypes = {'*.jpg;*.jpeg;*.png;*.tif;*.tiff',...
    'Image Files (*.jpg,*.jpeg,*.png,*.tif,*.tiff)'};
[customRefName, customRefPath] = uigetfile(imFileTypes, 'Pick custom reference');
customRef = imread(fullfile(customRefPath, customRefName));

[masterMapName, masterMapPath] = uigetfile(imFileTypes, 'Pick master map');
masterMap = imread(fullfile(masterMapPath, masterMapName));

% Crop master map

figure, imagesc(customRef), colormap gray
figure, imagesc(masterMap), colormap gray
h = imrect(gca, [1 1 512 512]); position = wait(h);

croppedMasterMap = imcrop(masterMap, position);

customRef = double(customRef); customRef = customRef./max(customRef(:));
croppedMasterMap = double(croppedMasterMap); croppedMasterMap = croppedMasterMap./max(croppedMasterMap(:));

dstr = datestr(now, 'mm_dd_yy_HH_MM_SS');

% On PC
saveRoot = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\ConeTransfer';

coneTransferFolder = ['coneTransfer_' dstr];

savePath = fullfile(saveRoot, coneTransferFolder);

if ~isdir(savePath)
    mkdir(savePath)
else
end


% On ICANDI

%savePath = '';
imwrite(customRef, fullfile(savePath, ['customRef_' dstr '.tif']));
imwrite(croppedMasterMap, fullfile(savePath, ['croppedMasterMap_' dstr '.tif']));
save(fullfile(savePath, ['rect_' dstr '.mat']), 'position');

end