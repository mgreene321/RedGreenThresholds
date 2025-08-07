%Load in data
subjectId = '20256R';
data_matrix = exp_data.data_matrix;
loc1_red = data_matrix(:,1) == 1 & data_matrix(:,3) == 1;

vidsOfInterest = find(loc1_red);

% Gather the frames from each video that are likely to contain cross
framesToLookAt = 13:19;

for i = 1:length(vidsOfInterest)
    vidObj = VideoReader(['20114R_' num2str(vidsOfInterest(i), '%03.f') '.avi']);
    vidFrames = read(vidObj);
    vidFrames = squeeze(vidFrames);
    vidFrames = vidFrames(:,:, framesToLookAt);
    allFrames(:,:,:,i) = vidFrames;
end

% For each frame of interest, locate the stimulus cross, and get the small
% region around the cross (lets say 100 x 100 pixels)

cross = zeros(11);
cross(6,:) = 1;
cross(:,6) = 1;

%Get cross locations
i = 1;
for v = 1:size(allFrames,4) %for each video   
    for f = 1:size(allFrames,3) % for each frame      
        frame = allFrames(:,:, f, v);   
        temp_frame = frame .* uint8(frame == 253); %red cross)        
        filteredFrame = filter2(cross, temp_frame);        
        [y x] = find(filteredFrame == max(filteredFrame(:)));        
        if length(x) > 1 || length(y) > 1
            cross_XY(i,:) = [nan nan f v];            
        else
            horizontalBarLocation = frame(y, x-5:x+5);
            if sum(horizontalBarLocation == 253) <  10                
                cross_XY(i,:) = [nan nan f v];
            else
                cross_XY(i,:) = [x y f v];
            end            
        end 
        i = i + 1;  
    end  
end

cross_XY(any(isnan(cross_XY),2),:) = [];

for v = min(cross_XY(:,4)):max(cross_XY(:,4))    
    for f = min(cross_XY(:,3)):max(cross_XY(:,3))        
        patch(:,:, f, v) = zeros(512, 512);                
        if ~isempty(cross_XY(cross_XY(:,3) == f & cross_XY(:,4) == v))            
            x = cross_XY(cross_XY(:,3) == f & cross_XY(:,4) == v, 1);
            y = cross_XY(cross_XY(:,3) == f & cross_XY(:,4) == v, 2);           
            frame = allFrames(:,:,f,v);          
            try
                patch(256-50:256+50, 256-50:256+50, f,v) = frame(y-50:y+50, x-50:x+50);
            catch
            end          
        else
        end        
    end    
end


%Compute mean image

%sum Frame
sumPatch = zeros(size(patch,1), size(patch,2));
nonZeroPatches = 0;
for v = min(cross_XY(:,4)):max(cross_XY(:,4))
    for f = min(cross_XY(:,3)):max(cross_XY(:,3))        
        if ~all(patch(:,:,f,v) == 0, "all")
            sumPatch = sumPatch + double(patch(:,:,f,v));
            nonZeroPatches = nonZeroPatches + 1;
        else
        end  
    end
end
meanPatch = sumPatch./nonZeroPatches;

%find how different each patch is from mean


%for each patch, find cone centers
cone_centers = zeros(size(patch));
for v = min(cross_XY(:,4)):max(cross_XY(:,4))
    
    for f = min(cross_XY(:,3)):max(cross_XY(:,3))
        
        if ~all(patch(:,:,f,v) == 0, "all")
            
            filtered_frame = filter2(Circle(3), patch(:,:,f,v).^2);
            filtered_frame = filtered_frame./max(filtered_frame(:));
            
            %FIND CONES
            
            cone_size = 5;
            xcorr_threshold = 0.6;
            
            [x_interest, y_interest] = img.find_cones(cone_size, filtered_frame, 'auto', xcorr_threshold,0);
            
            for i = 1:length(x_interest)
                cone_centers(y_interest(i), x_interest(i),f,v) = 1;
            end
            
        else
            
        end
    end
end

[x_interest, y_interest] = img.find_cones(cone_size, mean_patch, 'auto', xcorr_threshold,0);
mean_patch_cone_centers = zeros(512,512);
for i = 1:length(x_interest)
    mean_patch_cone_centers(y_interest(i), x_interest(i)) = 1;
end

filtered_mean_centers = imgaussfilt(mean_patch_cone_centers,3);
for i = 1:length(cross_XY)
    
    f = cross_XY(i,3); v = cross_XY(i,4);
    
    filtered_centers = imgaussfilt(cone_centers(:,:,f,v),3);
    
    similarity(i) = max(max(xcorr2(filtered_centers, filtered_mean_centers)));
    
end

%recompute mean_patch
patch_similarity = [cross_XY(:, [3 4]) similarity'];
patch_similarity = sortrows(patch_similarity,3);


weightedPatch = zeros(size(patch,1), size(patch,2));
for i = 1:length(patch_similarity)
    if similarity(i) > median(similarity)
   weightedPatch = weightedPatch + i.^2.*patch(:,:, patch_similarity(i,1), patch_similarity(i,2));
    else
    end
end
