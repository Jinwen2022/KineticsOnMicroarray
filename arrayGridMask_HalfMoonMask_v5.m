function [mask,topLeftCorner, spotCentroids,HalfMoonSpotsIDs, bkgMask] = arrayGridMask_HalfMoonMask_v5(imHalfMoon,imEq,pxSize,topLeftCorner,thresholdDNAIndices,badSpots,debugFlag,contrast)
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
%
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
radius = round(rowSpacing/3.9);
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
se = strel('disk',radius,6);
maskIni = imdilate(ceedMask,se);

diskSe = strel('disk',round(1.3*radius),6);
outMask = imdilate(ceedMask,diskSe);
rectSe = strel('rectangle',odd([round(rowSpacing*0.90) size(diskSe.Neighborhood,1)+3]));
bkgMaskIni = imdilate(ceedMask,rectSe);
bkgMaskIni = bkgMaskIni-outMask;

if debugFlag
    figure, imshow(labeloverlay(imEq,bkgMaskIni, 'Transparency', 0.3),contrast), hold on
    B = bwboundaries(maskIni,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Initial mask for computing binding intensities and selecting stronger binders than selected DNA operator')
end

[tmpFluoValues,tmpBkg] = computeSpotIntensities(imEq,maskIni,bkgMaskIni);%computed temFluoValues already is reducted by local background: tmpBkg
tmpFluoValues(badSpots) = NaN;
tmpBkg(badSpots) = NaN;
%Mask the spots that are most likely to have HalfMoon (unhomogenous
%binding) pattern, we selected spots that gave stronger binding than O1
radius_HalfMoon = round(2*radius);
se_HalfMoon = strel('disk',radius_HalfMoon,6);
mask_HalfMoon_Loose = imdilate(ceedMask,se_HalfMoon);
susHalfMoonSpotIDs = find(tmpFluoValues> mean(tmpFluoValues(thresholdDNAIndices)));
% Find the positions in mask_HalfMoon_Loose where the values match any of the values in susHalfMoonSpotIDs
indicesToRemove = ~ismember(mask_HalfMoon_Loose, susHalfMoonSpotIDs);
mask_HalfMoon_Loose(indicesToRemove) = 0;
temp_im =imEq;
temp_im(~mask_HalfMoon_Loose) = nanmedian(tmpBkg(:));
contrast = [nanmedian(tmpBkg(:)), nanmedian(tmpBkg(:))+ 1.25*nanmax(tmpFluoValues(:))];
temp_im2 = mat2gray(temp_im, contrast);
if debugFlag
    figure, imshow(temp_im2);
    title('Only stronger binders than selected DNA operator is showed')
end


% Step 2: 
% Finding the exact spot circle shape in the masked im with the strong binders mask
% generated just above
[centers, radii, ~] = imfindcircles(temp_im2, [round(0.9*radius) round(1.9*radius)], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);
%Remove indices in ceedMask that are not detected as strong binders
medianRadius = round(median(radii));  
se = strel('disk', medianRadius, 6);
circlesMask = false(size(ceedMask));
circleIndices = sub2ind(size(circlesMask), round(centers(:,2)), round(centers(:,1)));
circlesMask(circleIndices) = 1;
circlesMaskDilated = imdilate(circlesMask, se);%Mask out the regions that are detected as DNA spot for strong binders
% Initialize a mask to store the projected values
HalfMoonSpotsIDs = maskIni(circlesMask);
mask_HalfMoon_FitEdge = mask_HalfMoon_Loose;
mask_HalfMoon_FitEdge(~circlesMaskDilated) = 0;
if debugFlag
    figure, imshow(imEq,contrast), hold on
    B = bwboundaries(mask_HalfMoon_FitEdge,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Exact DNA spot for stronger binders  after circule detection with specific radius are masked out');
end


% Step 3: 
% Use detected circles centers for strong binding spots as anchors to
% transform  maskIni and bkgMaskIni for masking DNA microarray more
% accurately 
% Extract centroids and spot IDs from mask_HalfMoon_FitEdge. Pixel values
% will be used on step 4.
strongSpotProps = regionprops(mask_HalfMoon_FitEdge,imHalfMoon,'Centroid','Area','PixelValues','PixelList');
strongSpotLabels = find([strongSpotProps.Area]>0);
strongSpotProps = strongSpotProps(strongSpotLabels);
strongSpotCentroids = cat(1, strongSpotProps.Centroid);
% Extract corresponding centroids from maskIni (matching strong binder IDs)
ceedMaskProps = regionprops(ceedMask, 'Centroid');
ceedMaskCentroids = cat(1, ceedMaskProps.Centroid);
% Get the centroids of the corresponding strong binder IDs in maskIni
matchingceedMaskCentroids = ceedMaskCentroids(strongSpotLabels, :);
% Estimate the geometric transformation (affine) from the centroids
transformation = fitgeotrans(matchingceedMaskCentroids, strongSpotCentroids, 'affine');
% Apply the transformation to ceedMask to align it with im
outputRef = imref2d(size(imHalfMoon));  % Reference the size of image
ceedMaskTransformed = imwarp(ceedMask, transformation, 'OutputView', outputRef, 'InterpolationMethod', 'nearest');
% Generate the corrected mask
correctedMask = imdilate(ceedMaskTransformed, se);
if debugFlag
    figure, imshow(labeloverlay(imHalfMoon,maskIni, 'Transparency', 0.3),contrast), hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title({'Corrected mask(colorful labeloverlay) using transformation',...
        'based on detected strong binders positions',...
        'VS old mask (red boundries)'});
end
outMask = imdilate(ceedMaskTransformed,diskSe);
bkgMaskCorrected = imdilate(ceedMaskTransformed,rectSe);
bkgMaskCorrected = bkgMaskCorrected-outMask;
if debugFlag
    figure, imshow(labeloverlay(imHalfMoon,bkgMaskCorrected),contrast),hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('transformed maskIni with strongspots as anchors with corresponding corrected local bkg mask ');
end


% Step 4:
% Detect and mask top 20% binding are for each strong binding spot
% Initialize a new mask to store the top 20% regions on strong binding
% spots
maskTop20 = zeros(size(mask_HalfMoon_FitEdge));
% Loop through each unique label (each ROI)
for i = 1:length(strongSpotLabels)
    % Get the current label (ROI)
    label = strongSpotLabels(i);
    p = strongSpotProps(i);
    % Extract pixel intensities within this ROI from the image
    % Determine the threshold for the top 20% of intensities
    threshold = prctile(p.PixelValues, 80);  % Top 20% threshold (80th percentile)
    % Keep only the pixels that are above the threshold within this ROI
    top20Mask = p.PixelList(p.PixelValues >= threshold,:);
    % Retain the label in the new mask for the top 20% pixels
    maskTop20(sub2ind(size(maskTop20),top20Mask(:,2),top20Mask(:,1))) = label;
end
% use morphological operations for smoothing
smoothedTop20Mask = maskTop20;
se = strel('disk', round(medianRadius/5.3));  
smoothedTop20Mask = imdilate(smoothedTop20Mask, se);  % Dilate the mask to smooth
smoothedTop20Mask = imerode(smoothedTop20Mask, se);   % Erode to return to original size after dilation
% It's (almost) the same as smoothedTop20Mask=imclose(smoothedTop20Mask,se)
if debugFlag
    figure, imshow(labeloverlay(imadjust(imHalfMoon),maskTop20),contrast), hold on
    B = bwboundaries(correctedMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Smoothed top 20% pixels mask for only stronger binders than O1')
end


% Step 5:
% Compute centroids and labels for both strong and weak spots 
spotProps = regionprops(correctedMask, 'Centroid');
spotCentroids = round(cat(1, spotProps.Centroid));
% Find weak binding spots that are not in the strong spots
weakSpotLabels = 1:numel(spotProps);
weakSpotLabels(strongSpotLabels)=[];
% Compute distances between all weak and strong centroids
distancesMatrix = pdist2(spotCentroids(weakSpotLabels, :), strongSpotCentroids);
[~,nearestIdxs]=min(distancesMatrix,[],2);
% Initialize final mask
finalMask = smoothedTop20Mask;
% Crop in advance all strong spot masks
strongSpotMasks = zeros(2*medianRadius+1,2*medianRadius+1,length(strongSpotLabels));
for i = 1:length(strongSpotLabels)
    c = strongSpotCentroids(i,:);
    strongSpotMasks(:,:,i) = imcrop(smoothedTop20Mask, [c(1)-medianRadius,c(2)-medianRadius, 2*medianRadius,2*medianRadius]);
end
% Loop over each weak spot to copy the nearest strong spot mask
for i = 1:length(weakSpotLabels)
    weakLabel = weakSpotLabels(i);
    weakSpotCentroid = spotCentroids(weakLabel,:);
    % Find the nearest strong binding spot using precomputed distances matrix
    % Crop the mask for the nearest strong spot (top 20% mask)
    nearestStrongSpotHalfMoonMask = strongSpotMasks(:,:,nearestIdxs(i));
    transformedStrongSpotMask = nearestStrongSpotHalfMoonMask;
    transformedStrongSpotMask(nearestStrongSpotHalfMoonMask>0) = weakLabel;
    % Define weak spot dimensions
    weakXsCols = max(1,weakSpotCentroid(1)-medianRadius):min(size(finalMask, 2),weakSpotCentroid(1)+medianRadius);
    weakYsRows = max(1,weakSpotCentroid(2)-medianRadius):min(size(finalMask, 1),weakSpotCentroid(2)+medianRadius);
    % Resize the strong spot mask to fit within the weak spot dimensions
    if any([length(weakYsRows), length(weakXsCols)]<2*medianRadius+1)
        transformedStrongSpotMask = imresize(transformedStrongSpotMask, [length(weakYsRows), length(weakXsCols)], 'nearest');
    end
    % Assign the resized strong spot mask to the finalMask at the weak spot location
    finalMask(weakYsRows, weakXsCols) = transformedStrongSpotMask;
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