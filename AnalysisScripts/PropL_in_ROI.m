function PropL_in_ROI
%coneClassPath = 'C:\Programs\AOVIS_exp\+Max_exp\Cone_classifications';
%addpath(genpath(coneClassPath));

subjectId = '10001R';

coneSpacingPixels = 11;
    coneDiamPixels = coneSpacingPixels;
    coneAperturePixels = coneDiamPixels.*.5;
    gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
    windowSize = 11;
    coneAperture = fspecial('gaussian', windowSize, gaussianSigma);

    

%coneApertures = lightCaptureModel(subjectId);
%transferDataFolder = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\1_9_2024\ConeTransfer';
%coneApertures = dailyConeApertures(subjectId, transferDataFolder);

load('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\12_21_2023\TransferFromOCT\coneApertures.mat');
%coneApertures(:,:,1) = max(coneApertures(:)).*coneApertures(:,:,1)./max(coneApertures(:,:,:),[], 'all');
%coneApertures(:,:,2) = max(coneApertures(:)).*coneApertures(:,:,2)./max(coneApertures(:,:,:),[], 'all');
%coneApertures(:,:,3) = max(coneApertures(:)).*coneApertures(:,:,3)./max(coneApertures(:,:,:),[], 'all');

figure, imagesc(coneApertures./max(coneApertures(:))); hold on; axis equal;

xlim([0 size(coneApertures,2)]); ylim([0 size(coneApertures,1)]);

h = imrect(gca,[size(coneApertures,2)/2-10 size(coneApertures,1)/2-10, 21 21]);

setResizable(h,false);
addNewPositionCallback(h, @(pos) getPropL(pos, coneApertures));

function getPropL(pos, coneApertures)
pos = round(pos);
canvas = zeros(size(coneApertures,1), size(coneApertures,2));
canvas(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)) = 1;

mosaicUnderGrid = canvas .* coneApertures;

L = sum(mosaicUnderGrid(:,:,1), 'all', 'omitnan');
M = sum(mosaicUnderGrid(:,:,2), 'all', 'omitnan');
S = sum(mosaicUnderGrid(:,:,3), 'all', 'omitnan');

propL = L./(L+M+S);
%title([num2str(100*propL,'%0.2f') '% L, (' num2str(pos(1)) ', ' num2str(pos(2)) ') ' 'L+M+S= ' num2str(L+M+S)])
title(['L: ' num2str(L, '%0.2f') ', M: ' num2str(M, '%0.2f') ', S: ' num2str(S, '%0.2f')]);