function [mask,topLeftCorner, spotCentroids,HalfMoonSpotsIDs, bkgMask] = arrayGridMask_HalfMoonMask_v5(imHalfMoon,imEq,pxSize,topLeftCorner,thresholdDNAIndices,badSpotsMask,debugFlag,contrast)
% % The arrayGridMask_HalfMoonMask_v4 function creates a labeled mask for microarray spots, 
% focusing on masking the top 20% of pixels for strong binders and applying similar masks to weak binding spots. 
% It optionally generates a labeled mask for local background subtraction based on detected spot locations.
% 
% Inputs:
% 
% imHalfMoon: Preprocessed fluorescence image of the microarray.
% pxSize: Pixel size in microns.
% topLeftCorner: (Optional) Coordinates of the top-left spot centroid; if not provided, it can be set via GUI.
% thresholdDNAIndices: Indices of selected DNA operators that serve as a threshold to define strong binders.
% badSpotsMask: Array index for discarding unwanted Microarray REGIONS
% debugFlag: (Optional) Boolean flag to enable or disable debugging visualizations.
% contrast: (Optional) Intensity value cut-offs for image plotting, defaulting to [0, 660].


% Outputs:
% 
% mask: A refined mask image with segmented and labeled microarray spots.
% topLeftCorner: (x, y) coordinates of the top-left microarray spot centroid.
% spotCentroids: Centroid coordinates of all microarray spots.
% HalfMoonSpotsIDs: Indices of detected strong binding DNA spots.
% bkgMask: A labeled mask for background pixels around the microarray spots.
% Elf lab Jinwen Yuan 2024-09-24

% Step1:
colSpacing = 0.127*1000/pxSize;
rowSpacing = 0.073323*1000/pxSize;
radius = round(0.9*10*pxSize);
if size(imHalfMoon,1)>rowSpacing*265 && size(imHalfMoon,2)>colSpacing*84
    NROWS = 266;
    NCOLS = 85;
elseif size(imHalfMoon,1)>rowSpacing*95 && size(imHalfMoon,2)>colSpacing*81
    NROWS = 96;
    NCOLS = 82;
else
    error('The image size is too small')
end

if nargin<3 || isempty(topLeftCorner)
    figure, imshow(imHalfMoon(1:round(10*rowSpacing),1:round(10*colSpacing)),[0,660]);
    [x,y] = ginput(1);
    topLeftCorner = [x y];
end
if nargin<7
    debugFlag = false;
end
if nargin<8
    contrast = [0 660];
end

spotCentroidsX = topLeftCorner(1):colSpacing:colSpacing*NCOLS;
spotCentroidsY = topLeftCorner(2):rowSpacing:rowSpacing*NROWS;
spotCentroids = zeros(NROWS,NCOLS*2,2);

ceedMask = zeros(size(imHalfMoon));
[XX, YY] = meshgrid(spotCentroidsX, spotCentroidsY);
spotCentroids(:,1:2:end,1) = XX;
spotCentroids(:,1:2:end,2) = YY;
rcenters = round([XX(:) YY(:)]);
[XX, YY] = meshgrid(spotCentroidsX+colSpacing/2, spotCentroidsY+rowSpacing/2);
spotCentroids(:,2:2:end,1) = XX;
spotCentroids(:,2:2:end,2) = YY;
rcenters = [rcenters; round([XX(:) YY(:)])];
[~,ind] = sort(rcenters(:,1));
rcenters = rcenters(ind,:);
ceedMask(sub2ind(size(ceedMask),rcenters(:,2),rcenters(:,1))) = 1:numel(spotCentroids)/2;
se = strel('disk',radius,8);
maskIni = imdilate(ceedMask,se);

if debugFlag
    figure, imshow(imEq,contrast), hold on
    B = bwboundaries(maskIni,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Initial mask for computing binding intensities and selecting stronger binders than selected DNA operator')
end

diskSe = strel('disk',round(1.35*10*pxSize),8);
outMask = imdilate(ceedMask,diskSe);
rectSe = strel('rectangle',odd([round(rowSpacing*0.95) size(diskSe.Neighborhood,1)+3]));
bkgMaskIni = imdilate(ceedMask,rectSe);
bkgMaskIni = bkgMaskIni-outMask;
if debugFlag
    figure, imshow(labeloverlay(imadjust(imEq),bkgMaskIni)),axis image
end

[tmpFluoValues,tmpBkg] = computeSpotIntensities(imEq,maskIni,bkgMaskIni);%computed temFluoValues already is reducted by local background: tmpBkg

%Mask the spots that are most likely to have HalfMoon (unhomogenous
%binding) pattern, we selected spots that gave stronger binding than O1
radius_HalfMoon = round(0.036662*1000/pxSize);
se_HalfMoon = strel('disk',radius_HalfMoon,8);
mask_HalfMoon_Loose = imdilate(ceedMask,se_HalfMoon);
mask_HalfMoon_Loose(badSpotsMask) = 0;
susHalfMoonSpotIDs = find(tmpFluoValues> mean(tmpFluoValues(thresholdDNAIndices)));
% Find the positions in mask_HalfMoon_Loose where the values match any of the values in susHalfMoonSpotIDs
indicesToRemove = ~ismember(mask_HalfMoon_Loose, susHalfMoonSpotIDs);
indicesToKeep = ismember(mask_HalfMoon_Loose, susHalfMoonSpotIDs);
mask_HalfMoon_Loose(indicesToRemove) = 0;
tempMask = mask_HalfMoon_Loose;
tempMask(indicesToKeep) = 1;
temp_im =imEq;
temp_im(~tempMask) =median(median(tmpBkg));
if debugFlag
    figure, imshow(imEq,contrast), hold on
    B = bwboundaries(tempMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Only stronger binders than selected DNA operator is masked')
end


% Step 2: 
% Finding the exact spot circle shape in the masked im with the strong binders mask
% generated just above
temp_im2 = imadjust(temp_im);
[centers, radii, ~] = imfindcircles(temp_im2, [round(radius) round(0.9*radius_HalfMoon)], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);
%Remove indices in ceedMask that are not detected as strong binders
medianRadius = round(median(radii));  
se = strel('disk', medianRadius, 8);
circlesMask = logical(zeros(size(ceedMask)));
circleIndices = sub2ind(size(circlesMask), round(centers(:,2)), round(centers(:,1)));
circlesMask(circleIndices) = 1;
circlesMaskDilated = imdilate(circlesMask, se);%Mask out the regions that are detected as DNA spot for strong binders
% Initialize a mask to store the projected values
ceedMaskFiltered = zeros(size(ceedMask));
% project values from mask onto positions labelled by circlesMask
ceedMaskFiltered(circlesMask) = maskIni(circlesMask);  
HalfMoonSpotsIDs = maskIni(circlesMask);
mask_HalfMoon_FitEdge=mask_HalfMoon_Loose;
mask_HalfMoon_FitEdge(~circlesMaskDilated) = 0;
if debugFlag
    figure, imshow(imEq,[0,1000]), hold on
    B = bwboundaries(mask_HalfMoon_FitEdge,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Exact DNA spot for stronger binders  after circule detection with specific radius are masked out');
end


% Step 3: 
% Use detected circles centers for strong binding spots as anchors to
% transform  maskIni and bkgMaskIni for maksing DNA microarray more
% accurately 
% Extract centroids and spot IDs from mask_HalfMoon_FitEdge
propsHalfMoon = regionprops(mask_HalfMoon_FitEdge, 'Centroid');
strongBinderCentroids = cat(1, propsHalfMoon.Centroid);
validStrongBinderIdx = ~any(isnan(strongBinderCentroids), 2);  % Identify valid centroids
strongBinderCentroids = strongBinderCentroids(validStrongBinderIdx, :);  % Keep only valid centroids
strongBinderIDs = unique(mask_HalfMoon_FitEdge(mask_HalfMoon_FitEdge > 0)); % Only non-zero strong binder spots
% Extract corresponding centroids from maskIni (matching strong binder IDs)
ceedMaskProps = regionprops(ceedMask, 'Centroid');
ceedMaskCentroids = cat(1, ceedMaskProps.Centroid);
% Get the centroids of the corresponding strong binder IDs in maskIni
matchingceedMaskCentroids = ceedMaskCentroids(strongBinderIDs, :);
% Estimate the geometric transformation (affine) from the centroids
transformation = fitgeotrans(matchingceedMaskCentroids, strongBinderCentroids, 'affine');
%Apply the transformation to ceedMask to align it with im
outputRef = imref2d(size(imHalfMoon));  % Reference the size of image
ceedMaskTransformed = imwarp(ceedMask, transformation, 'OutputView', outputRef, 'InterpolationMethod', 'nearest');
% Generate the corrected mask
correctedMask = imdilate(ceedMaskTransformed, se);
% correctedMask(mask_HalfMoon_FitEdge > 0) = mask_HalfMoon_FitEdge(mask_HalfMoon_FitEdge > 0); % Replace strong binder areas with precise mask
if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon),maskIni, 'Transparency', 0.3),contrast), hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title({['Corrected mask(colorful labeloverlay) using transformation'],...
        ['based on detected strong binders positions'],...
        ['VS old mask (red boundries)']});
end
diskSe = strel('disk',round(1.35*10*pxSize),8);
outMask = imdilate(ceedMaskTransformed,diskSe);
rectSe = strel('rectangle',odd([round(rowSpacing*0.95) size(diskSe.Neighborhood,1)+3]));
bkgMaskCorrected = imdilate(ceedMaskTransformed,rectSe);
bkgMaskCorrected = bkgMaskCorrected-outMask;
if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon),bkgMaskCorrected)),hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('transformed maskIni with strongspots as anchors with corresponding corrected local bkg mask ');
end



% Step 4:
% Detect and mask top 20% binding are for each strong binding spot

% Initialize a new mask to store the top 20% regions on strong binding
% spots
maskTop20 = zeros(size(mask_HalfMoon_FitEdge));
mask = correctedMask;
% Get all unique labels in the mask (excluding zero)
strongSpotLabels = unique(mask_HalfMoon_FitEdge(mask_HalfMoon_FitEdge > 0));
% Loop through each unique label (each ROI)
for i = 1:length(strongSpotLabels)
    % Get the current label (ROI)
    label = strongSpotLabels(i);
    % Extract the region corresponding to this label from the mask
    roiMask = (mask_HalfMoon_FitEdge == label);
    spotMask_mask = (mask == label);
    % Extract pixel intensities within this ROI from the image
    roiIntensities = imHalfMoon(roiMask);
    % Determine the threshold for the top 20% of intensities
    threshold = prctile(roiIntensities, 80);  % Top 20% threshold (80th percentile)
    % Keep only the pixels that are above the threshold within this ROI
    top20Mask = roiMask & (imHalfMoon >= threshold);
    % Retain the label in the new mask for the top 20% pixels
    maskTop20(top20Mask) = label;
    % Remove label's region in original mask, replace with new top 20% mask
    mask(spotMask_mask) = 0;
    mask(top20Mask) = label;
end
% use morphological operations for smoothing
smoothedTop20Mask = maskTop20;
se = strel('disk', 3);  
smoothedTop20Mask = imdilate(smoothedTop20Mask, se);  % Dilate the mask to smooth
smoothedTop20Mask = imerode(smoothedTop20Mask, se);   % Erode to return to original size after dilation
if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon),smoothedTop20Mask),contrast), hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Smoothed top 20% pixels mask for only stronger binders than O1')
end




% Step 5:
% Compute centroids and labels for both strong and weak spots 
strongSpotProps = regionprops(mask_HalfMoon_FitEdge>0, "Centroid");
strongSpotsCentroids = round(cat(1, strongSpotProps.Centroid));
strongSpotsLabels = mask_HalfMoon_FitEdge(sub2ind(size(mask_HalfMoon_FitEdge),strongSpotsCentroids(:,2),strongSpotsCentroids(:,1)));
spotProps = regionprops(correctedMask, 'Centroid');
spotCentroids = round(cat(1, spotProps.Centroid));
spotIDs = correctedMask(sub2ind(size(correctedMask),spotCentroids(:,2),spotCentroids(:,1)));
% Find weak binding spots that are not in the strong spots
[weakSpotIDs,indexWeakSpotIDs] = setdiff(spotIDs, strongSpotsLabels); %spotIDs (indexSpotIDs) = weakSpotIDs
% Compute distances between all weak and strong centroids
distancesMatrix = pdist2(spotCentroids(indexWeakSpotIDs, :), strongSpotsCentroids);%Note strongSpotLabels is not correctly posisionted related to strongSpotCentriods
% Initialize final mask
finalMask = smoothedTop20Mask;
% Loop over each weak spot to copy the nearest strong spot mask
for i = 1:length(weakSpotIDs)
    weakLabel = weakSpotIDs(i);
    weakSpotCentroid = spotCentroids(indexWeakSpotIDs(i),:);
    % Find the nearest strong binding spot using precomputed distances matrix
    [~, nearestIdx] = min(distancesMatrix(i, :));
    nearestStrongspotCentroid = strongSpotsCentroids (nearestIdx,:);
    % Crop the mask for the nearest strong spot (top 20% mask)
    nearestStrongSpotHalfMoonMask = imcrop(smoothedTop20Mask, [nearestStrongspotCentroid(1)-medianRadius,nearestStrongspotCentroid(2)-medianRadius, 2*medianRadius,2*medianRadius]);
    transformedStrongSpotMask = double(logical(nearestStrongSpotHalfMoonMask))*weakLabel;
    % Define weak spot dimensions
    weakXsCols = weakSpotCentroid(1)-medianRadius:weakSpotCentroid(1)+medianRadius;
    weakYsRows = weakSpotCentroid(2)-medianRadius:weakSpotCentroid(2)+medianRadius;
    % Adjust rows and cols to fit within image bounds
    validWeakXsCols = weakXsCols(weakXsCols > 0 & weakXsCols <= size(finalMask, 2));
    validWeakYsRows = weakYsRows(weakYsRows > 0 & weakYsRows <= size(finalMask, 1));
    % Resize the strong spot mask to fit within the weak spot dimensions
    transformedStrongSpotMask = imresize(transformedStrongSpotMask, [length(validWeakYsRows), length(validWeakXsCols)], 'nearest');
    % Assign the resized strong spot mask to the finalMask at the weak spot location
    finalMask(validWeakYsRows, validWeakXsCols) = transformedStrongSpotMask;
end


mask = finalMask;
bkgMask = bkgMaskCorrected;

if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon), finalMask), contrast), hold on
    B = bwboundaries(smoothedTop20Mask, 'noholes');
    visboundaries(B, 'EnhanceVisibility', 0, 'LineWidth', 1);
    title('Final Mask with Top 20% for Strong Spots and Copied nearest strong spots for Weak Spots');
end

if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon),bkgMask),contrast), hold on
    B = bwboundaries(maskTop20,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    B = bwboundaries(mask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Masking out top 20% pixels among every DNA spot with their local background');
end



end