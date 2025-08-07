%Load in the coneApertures matrix, which is an array of red, green, and
%blue 2D Gaussians that represent the positions and aperture functions of
%L, M, and S-cones

tic;
coneApertures = coneApertures./max(coneApertures(:));
L = coneApertures(:,:,1);
M = coneApertures(:,:,2);

PPD = 560;
sideLength = 2.25; % arcmin
gridWidth = 2;
gridHeight = 2;
sideLength_pix = PPD .* sideLength ./ 60;
avgKernel = ones(round(sideLength_pix));
avgKernel = avgKernel ./ numel(avgKernel);

L_avg = conv2(L, avgKernel, 'same');
M_avg = conv2(M, avgKernel, 'same');

propL_avg = L_avg./(L_avg + M_avg);

[rows,cols] = find(sum(coneApertures,3));

minRow = round(min(rows) + 2.*sideLength_pix);
maxRow = round(max(rows) - 2.*sideLength_pix);
minCol = round(min(cols) + 2.*sideLength_pix);
maxCol = round(max(cols) - 2.*sideLength_pix);

%topLeftCellCenters = [topLeftCellCenters_Y(:) topLeftCellCenters_X(:)];

%For every grid position, get L/(L+M) in each cell and put it in matrix

% | x1 x2 |
% | x3 x4 |

propL_in_cells = [];
topLeftCellCenters = [];

LM_sum = L_avg + M_avg;
LM_sum = LM_sum(LM_sum > 0);

min_LM = median(LM_sum(:)) - 1.*std(LM_sum,1, 'all');
max_LM = median(LM_sum(:)) + 2.*std(LM_sum,1, 'all');


%Get every possible gridIdxation for the center of the top left cell center

[X Y] = meshgrid(minCol:maxCol, minRow:maxRow);
topLeftCellCenters = [X(:) Y(:)]; % scatter(topLeftCellCenters(i,1), topLeftCellCenters(i,2)) will plot this gridIdxation on coneApertures image
propL_in_cells = nan(size(topLeftCellCenters,1), gridWidth.*gridHeight);
removeIdx = zeros(size(topLeftCellCenters,1), 1);

parfor gridIdx = 1:size(topLeftCellCenters,1)
    [gridx gridy] =  meshgrid(0:gridWidth-1, 0:gridHeight-1);
    
    topLeftCellCenter = topLeftCellCenters(gridIdx,:);
    gridCellCenters_Y = topLeftCellCenter(2) + round(sideLength_pix).*gridy;
    gridCellCenters_X = topLeftCellCenter(1) + round(sideLength_pix).*gridx;
    propL_temp =  propL_avg(unique(gridCellCenters_Y(:)), unique(gridCellCenters_X(:)));
    propL_temp = transpose(propL_temp(:));
    
    propL_in_cells(gridIdx,:) = propL_temp;
    
    L_temp = L_avg(unique(gridCellCenters_Y(:)), unique(gridCellCenters_X(:)));
    M_temp = M_avg(unique(gridCellCenters_Y(:)), unique(gridCellCenters_X(:)));
    LM_temp = L_temp + M_temp;
    
    
    if any(any(isnan(propL_temp))) || ...
            any(any(LM_temp < min_LM)) || ...
            any(any(LM_temp > max_LM))
        removeIdx(gridIdx) = true;
    else
    end
end
toc;

propL_in_cells(find(removeIdx),:) = [];
topLeftCellCenters(find(removeIdx),:) = [];

for gridIdx = 1:size(topLeftCellCenters,1)
    yi = transpose(propL_in_cells(gridIdx,:));
    yi = sort(yi);
    xi(:,2) = transpose(1:4);
    xi(:,1) = 1;
    
    b = xi \ yi;
    yi_pred = transpose((1:4).*(b(2)) + b(1));
    
    Rsq = 1 - sum((yi - yi_pred).^2)./sum((yi - mean(yi).^2));
    gridStats(gridIdx, 1) = b(2);
    gridStats(gridIdx,2) = Rsq;
    
end

gridStatsSum = gridStats(:,1) + gridStats(:,2);
gridStatsSum(:,2) = 1:length(gridStatsSum);
gridStatsSum = sortrows(gridStatsSum,1, 'descend');

gridStatsSum(any(isnan(gridStatsSum),2),:) = [];

figure, imagesc(coneApertures); hold on
[gridx gridy] =  meshgrid(0:gridWidth-1, 0:gridHeight-1);
topLeftCellCenter = topLeftCellCenters(gridStatsSum(1,2),:);
gridCellCenters_Y = topLeftCellCenter(2) + round(sideLength_pix).*gridy;
gridCellCenters_X = topLeftCellCenter(1) + round(sideLength_pix).*gridx;

rectPos_Y = gridCellCenters_Y - round(sideLength_pix./2);
rectPos_X = gridCellCenters_X - round(sideLength_pix./2);

for i = 1:gridWidth
    for j = 1:gridHeight
        pos = [rectPos_X(j,i) rectPos_Y(j,i) round(sideLength_pix) round(sideLength_pix)];
        drawrectangle('Position', pos, 'LineWidth', 1, 'Color', 'w');
    end
end