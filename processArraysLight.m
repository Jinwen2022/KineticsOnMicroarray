function [fluoValues, mask] = processArraysLight(dirFluo,dirOffset,parameters)
% Performs analysis of specific operator binding experiments using
% microarrays without saving preprocessed images. Handles either Cy3 or Cy5
% images.
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
%
% Output:
%   fluoValues - 3D (or 2D for Cy5) array, where (f,i,j) is the intensity
%                of microarray spot (i,j) at acquisition frame f.
%   mask - labelled image of microarray spots in the preprocessed images.
%
% Spartak Zikrin, Elf lab, 2024-05-17.
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
imSize = parameters.roi([4 3])+1;
slide = zeros(imSize);
[mask, ~, ~, bkgMask] = arrayGridMask(slide,parameters.pxSize,parameters.topLeftCorner);
%% Preprocess Cy3 images and extract spot intensities
tmpValuesCy3 = cell(1,numel(dirFluo));
for i=1:numel(dirFluo)
    tmpValuesCy3{i} = computeSpotIntensitiesFromRawData(dirFluo{i},parameters.overlap,...
        parameters.angularDisplacementTile,parameters.rangePositions,parameters.posColRange,...
        parameters.weights,parameters.angularDisplacement,parameters.roi,offset,mask,bkgMask);
end
fluoValues = squeeze(vertcat(tmpValuesCy3{:}));