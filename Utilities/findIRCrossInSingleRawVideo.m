function IR_cross_XY = findIRCrossInSingleRawVideo(vidFileName, frameNums)

% subjectId = '10001R';
% expDate = '1_25_2024';
% expTime = '15_27_9';
% 
% RedGreenThresholdsPath = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\';
% goodDeliveries = importdata(fullfile(RedGreenThresholdsPath, subjectId, expDate, expTime, 'goodDeliveries.mat'));
% 
% vidPath = 'Z:\Local_Share\Imaging_Data\AOVIS\10001\Right\1_25_2024_15_27_9';
% vidFile =  ['10001R_' num2str(15,'%.3d'), '.avi'];

vidObj = VideoReader(vidFileName);


for i = 1:numel(frameNums)

    frame = read(vidObj, frameNums(i));
try
    IR_cross_XY(i,:) = findIRCross(frame);
catch
    IR_cross_XY(i,:) = [nan nan];
end

end

function IR_cross_XY = findIRCross(frame)

% Make stimulus cross
cross_diam = 11;
cross = zeros(cross_diam);
cross(6,:) = 1;
cross(:,6) = 1;

temp_IR_frame = frame .* uint8(frame == 255);

filtered_IR_frame = filter2(cross, temp_IR_frame);

[y_IR, x_IR] = find(filtered_IR_frame == max(filtered_IR_frame(:)));

% Generate black image
blank_im = zeros([size(frame) 3]);

if crossMightBeThere(x_IR, y_IR)
    x_IR = x_IR(1); y_IR = y_IR(1);
    if x_IR-floor(cross_diam/2) > 0 && y_IR-floor(cross_diam/2) > 0
        % Draw IR cross
        blank_im(y_IR, x_IR-floor(cross_diam/2):x_IR+floor(cross_diam/2),1) = 1;
        blank_im(y_IR-floor(cross_diam/2):y_IR+floor(cross_diam/2), x_IR,1) = 1;
    else
        x_IR = nan; y_IR = nan;
    end
else
    x_IR = nan; y_IR =nan;
end

hypoIRCross = xor(blank_im(:,:,1), sum(blank_im,3)>1) .* blank_im(:,:,1);
avgFilter = ones(3,3);
IRCrossSearch = conv2(hypoIRCross, avgFilter, 'same') > 0;

if sum(hypoIRCross(:)) > 1 && crossIsThere(x_IR, y_IR, frame, hypoIRCross, IRCrossSearch, 255)
    IR_cross_XY = [x_IR y_IR];
else
  
end


% Auxiliary functions

function bool = crossMightBeThere(x_cross, y_cross)

% When filtering image for cross, two nearby peaks may be found

if numel(x_cross) == 1 && numel(y_cross) == 1
    bool = true;
elseif numel(x_cross) == 2 && numel(y_cross) == 2

    if abs(diff(x_cross)) <= 1 && abs(diff(y_cross)) <= 1
        bool = true;
    else
        bool = false;
    end
else
    bool = false;

end

function bool = crossIsThere(x_cross, y_cross, frame, hypoCross, crossSearch, crossVal)

% Confirms the presence of a cross if either the entire horizontal bar is
% found, the entire vertical bar is found, or if the two halves of the
% veritcal bar are found slightly offset from one another.

maskedFrame = uint8(crossSearch).*frame;

% conditions to find cross

searchRegionFilled = sum(maskedFrame(:) == crossVal) >= sum(hypoCross(:)) - 1;
horizontalBarFound = sum(maskedFrame(y_cross,:) == crossVal) >= sum(hypoCross(y_cross,:)) - 1;
verticalBarFound =  sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal, 'all') >= sum(hypoCross(:,x_cross)) - 1 &&...
    ~any(sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal, 2) > 1); %entire, possibly sheared vertical bar found


horizontalBarNotTooTall = sum(sum(maskedFrame(y_cross-1:y_cross+1, :) == crossVal,1) > 2) <= 4;
verticalBarNotTooWide = sum(sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal,2) > 2) <= 4;

if horizontalBarFound && horizontalBarNotTooTall
    bool = true;
elseif verticalBarFound && verticalBarNotTooWide
    bool = true;
else
    bool = false;
end