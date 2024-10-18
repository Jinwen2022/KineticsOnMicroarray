function [fluoValues, bkgValues, mask,bkgMask,HalfMoonSpotsIDs] = processArraysLight_HalfMoonMask_v5(dirFluo,dirOffset,parameters,thresholdDNAIndices, badSpots)
% Performs analysis of specific operator binding experiments using
% microarrays without saving preprocessed images. Handles either Cy3 or Cy5
% images. 

%This function version 3 compute the top 20% Cy3values on strong binders'
% DNA spot, which is defined stronger as give stronger Cy3values than
% selected DNA (thresholdDNAIndices) in equilibrium with averaged spot
% Intensity values. Local background reduction is done by using original
% bkgMask
%
% Input:
%   dirFluo - Path to directory (or cell array with paths to directories)
%            with fluorescence images. Or a directory
%   dirOffset - path to directory with an offset image for
%               background subtraction. If empty no background subtraction.
%   parameters - scructure array with fields:
%       overlap - number of tile overlapping pixels.
%       angularDisplacementTile - in degrees, rotate counter-clockwise
%                                 raw images (tiles).
%       angularDisplacement - in degrees, rotate counter-clockwise
%                             merged images.
%       rangePositions - array with position indices. If empty,
%                           all positions are used. Default [].
%       roi - cropping region to crop images after rotation.
%       pxSize - pixel size in mum.
%       topLeftCorner - [x y] coordinates of the centroid of the
%                       (imaginary) spot in the top left corner of
%                       preprocessed image.
%   thresholdDNAIndices - chosen DNA operator's indexes in DNA array for
%                       computing strong binder's fluoValues
%   badSpotsMask - badSpots logical array generated from manuelly labelled badSpotsMask, to discrad unwanted
%   DNA spots.
%
% Output:
%   fluoValues - 3D (or 2D for Cy5) array, where (f,i,j) is the intensity
%                of microarray spot (i,j) at acquisition frame f.
%   bkgValues  - 3D (or 2D for Cy5) array, where (f,i,j) is the average 
%                pixel intensity in the neighborhood of microarray spot 
%                (i,j) at acquisition frame f.
%   HalfMoonSpotsIDs - 
%   mask - labelled image of spots suspicious of Unhomogenous binding in the preprocessed images.
%   bkgMask - labelled image of spots suspicious of Unhomogenous binding's local background in the preprocessed images
%
% Spartak Zikrin, Elf lab, 2024-05-17.
% Jinwen Yuan, Modified, 2024-09-19.
%% Parse parameters
if ~isfield(parameters,'rangePositions')
    parameters.rangePositions = [];
end
if ~isfield(parameters,'weights')
    parameters.weights = [];
end
if ~isfield(parameters,'posColRange')
    parameters.posColRange = [];
end
if ~iscell(dirFluo)
    dirFluo = {dirFluo};
end
%% Preprocess offset image if it is available
if ~isempty(dirOffset)
    offset = stitchArrayImages(dirOffset,parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,1);
    offset = preprocessImage(offset,parameters.angularDisplacement,parameters.roi);
else
    offset = 0;
end
%% Compute spot mask and background mask
slideHalfMoon = stitchArrayImages(dirFluo{1},parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,33);
slideEq = stitchArrayImages(dirFluo{2},parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,1);
imHalfMoon = preprocessImage(slideHalfMoon,parameters.angularDisplacement,parameters.roi);
imEq = preprocessImage(slideEq,parameters.angularDisplacement,parameters.roi);
[mask,~, ~,HalfMoonSpotsIDs, bkgMask] = arrayGridMask_HalfMoonMask_v5(imHalfMoon,imEq,parameters.pxSize,parameters.topLeftCorner,thresholdDNAIndices,badSpots);
%% Preprocess Cy3 images and extract spot intensities
tmpFluoValues = cell(1,numel(dirFluo));
tmpBkgValues = cell(1,numel(dirFluo));
for i=1:numel(dirFluo)
    [tmpFluoValues{i}, tmpBkgValues{i}] = computeSpotIntensitiesFromRawData(dirFluo{i},parameters.overlap,...
        parameters.angularDisplacementTile,parameters.rangePositions,parameters.posColRange,...
        parameters.weights,parameters.angularDisplacement,parameters.roi,offset,mask,bkgMask);
end
fluoValues = squeeze(vertcat(tmpFluoValues{:}));
bkgValues = squeeze(vertcat(tmpBkgValues{:}));