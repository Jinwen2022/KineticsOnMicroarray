function [maskIni,topLeftCorner, bkgMask,HalfMoonSpotsIDs, mask_HalfMoon_FitEdge] = arrayGridMask_EqMask(im,pxSize,topLeftCorner,thresholdDNAIndices,debugFlag,contrast)
% Creates a labelled mask of strong binders among microarray spots and optionally a labelled mask for local
% background subtraction.
% The spot mask is only 20% of the real spot are in order to deal with
% the unhomogenous binding pattern for strong binders, with bkgMask
% generated based on circle detection on strong binders.
%
% Input:
%   im - preprocessed fluo image.
%   pxSize - pixel size in mum.
%   topLeftCorner - [x,y] coordinates of the top left spot centroid.
%                   Optional, if not provided or empty, the centroid is
%                   set using gui.  Default [].
%   thresholdDNAIndices - selected certain DNA operator indices in DNA
%                   array, as the threshold to define strong operators that most likely
%                   causing unhomogenous binding on DNA spot.
%   debugFlag - optional, if true plot image with spot outlines. Default
%               false.
%   contrast - a tuple with upper and lower pixel intensity value cut-offs
%              for image plotting. By default, [0 660].

%
% Output:
%   mask - an image of the same size as im with segmented and labelled
%          microarray spots;
%   topLeftCorner - (x,y)-coordinates of the top left microarray spot
%                   centroid;
%   spotCentroids - coordinates of all microarray spot centroids
%                   HalfMoonSpotsIDs - indices of all deteced strong binders (stronger than
%                   the input DNA operator that is defined by thresholdDNAIndices)
%   HalfMoonSpotsIDs - Indices of detected strong binding DNA spots in dna Sequences
%                   printed on Microarray
%   bkgMask - optional, a labelled mask of background pixels around
%             microarray spots


colSpacing = 0.127*1000/pxSize;
rowSpacing = 0.073323*1000/pxSize;
radius = round(0.9*10*pxSize);
if size(im,1)>rowSpacing*265 && size(im,2)>colSpacing*84
    NROWS = 266;
    NCOLS = 85;
elseif size(im,1)>rowSpacing*95 && size(im,2)>colSpacing*81
    NROWS = 96;
    NCOLS = 82;
else
    error('The image size is too small')
end

if nargin<3 || isempty(topLeftCorner)
    figure, imshow(im(1:round(10*rowSpacing),1:round(10*colSpacing)),[0,660]);
    [x,y] = ginput(1);
    topLeftCorner = [x y];
end
if nargin<5
    debugFlag = false;
end
if nargin<6
    contrast = [0 660];
end

spotCentroidsX = topLeftCorner(1):colSpacing:colSpacing*NCOLS;
spotCentroidsY = topLeftCorner(2):rowSpacing:rowSpacing*NROWS;
spotCentroids = zeros(NROWS,NCOLS*2,2);

ceedMask = zeros(size(im));
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
    figure, imshow(im,contrast), hold on
    B = bwboundaries(maskIni,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
end

diskSe = strel('disk',round(1.35*10*pxSize),8);
outMask = imdilate(ceedMask,diskSe);
rectSe = strel('rectangle',odd([round(rowSpacing*0.95) size(diskSe.Neighborhood,1)+3]));
bkgMaskIni = imdilate(ceedMask,rectSe);
bkgMaskIni = bkgMaskIni-outMask;
if debugFlag
    figure, imshow(labeloverlay(imadjust(im),bkgMaskIni)),axis image
end

[tmpFluoValues,tmpBkg] = computeSpotIntensities(im,maskIni,bkgMaskIni);%computed temFluoValues already is reducted by local background: tmpBkg

%Mask the spots that are most likely to have HalfMoon (unhomogenous
%binding) pattern, we selected spots that gave stronger binding than O1
radius_HalfMoon = round(0.036662*1000/pxSize);
se_HalfMoon = strel('disk',radius_HalfMoon,8);
mask_HalfMoon_Loose = imdilate(ceedMask,se_HalfMoon);
susHalfMoonSpotIDs = find(tmpFluoValues> mean(tmpFluoValues(thresholdDNAIndices)));
% Find the positions in mask_HalfMoon_Loose where the values match any of the values in susHalfMoonSpotIDs
indicesToRemove = ~ismember(mask_HalfMoon_Loose, susHalfMoonSpotIDs);
indicesToKeep = ismember(mask_HalfMoon_Loose, susHalfMoonSpotIDs);
mask_HalfMoon_Loose(indicesToRemove) = 0;
tempMask = mask_HalfMoon_Loose;
tempMask(indicesToKeep) = 1;
temp_im =im;
temp_im(~tempMask) =median(median(tmpBkg));
if debugFlag
    figure, imshow(temp_im,contrast), hold on
    B = bwboundaries(tempMask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Mask1 for masking out stronger binders than selected DNA operator')
end


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
    figure, imshow(im,[0,1000]), hold on
    B = bwboundaries(mask_HalfMoon_FitEdge,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Masking out exact DNA spot for stronger binders than O1');
end


% Initialize a new mask to store the top 20% regions
maskTop20 = zeros(size(mask_HalfMoon_FitEdge));
mask = maskIni;
% Get all unique labels in the mask (excluding zero)
uniqueLabels = unique(mask_HalfMoon_FitEdge);
uniqueLabels(uniqueLabels == 0) = [];  
% Loop through each unique label (each ROI)
for i = 1:length(uniqueLabels)
    % Get the current label (ROI)
    label = uniqueLabels(i);
    % Extract the region corresponding to this label from the mask
    roiMask = (mask_HalfMoon_FitEdge == label);
    spotMask_mask = (mask == label);
    % Extract pixel intensities within this ROI from the image
    roiIntensities = im(roiMask);
    % Determine the threshold for the top 20% of intensities
    threshold = prctile(roiIntensities, 80);  % Top 20% threshold (80th percentile)
    % Keep only the pixels that are above the threshold within this ROI
    top20Mask = roiMask & (im >= threshold);
    % Retain the label in the new mask for the top 20% pixels
    maskTop20(top20Mask) = label;
    % Remove label's region in original mask, replace with new top 20% mask
    mask(spotMask_mask) = 0;
    mask(top20Mask) = label;
end

bkgMask = bkgMaskIni;

if debugFlag
    figure, imshow(labeloverlay(imadjust(im),bkgMask),contrast), hold on
    B = bwboundaries(maskTop20,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    B = bwboundaries(mask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
    title('Masking out top 20% pixels among every exact DNA spot for stronger binders than O1');
end



end