function [ValuesCy5, mask, bkg] = processArraysLight_Cy5_temp(dirCy5,dirOffset,parameters,badSpots)
% Performs analysis of specific operator binding experiments using
% microarrays. Only Cy3 images are being analysed and preprocessed images
% are not saved.
%
% Input:
%   dirCy3 - cell array, each element is a path to directory with Cy3
%            fluorescence images.
%   dirOffset - path to directory with an offset Cy3 image, used for
%               background subtraction. If empty no background subtraction.
%   parameters - scructure array with fields:
%       overlap - number of tile overlapping pixels.
%       angularDisplacementTile - in degrees, rotate counter-clockwise
%                                 raw images (tiles).
%       angularDisplacement - in degrees, rotate counter-clockwise
%                             merged Cy3 images.
%       rangePositions - array with position indices in Cy3. If empty,
%                           all positions are used. Default [].
%       roi - cropping region in Cy3 images after rotation.
%       pxSize - pixel size in mum.
%       topLeftCorner - [x y] coordinates of the centroid of the
%                       (imaginary) spot in the top left corner of
%                       preprocessed image.
%   badSpots - mask of spots to be excluded for computing average
%              background. Optional, default [].
%
% Output:
%   valuesCy3 - cell 2D array, each element corresponds to a microarray
%               spot and contains a vector of spot intensities in Cy3
%               sorted by image acquisition time.
%   mask - boolean matrix of microarray spots in the preprocessed images.
%   bkg - array with median background intensity. Optional.
%
% Spartak Zikrin, Elf lab, 2021-09-12.
%% Parse parameters
if ~isfield(parameters,'rangePositions')
    parameters.rangePositions = [];
end
if nargin<4
    badSpots = [];
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
mask = arrayGridMask_temp(slide,parameters.pxSize,parameters.topLeftCorner);
if nargout==3
    computeBkg = true;
    bkgMask = false(imSize);
    rowSpacing = 0.073323*1000/parameters.pxSize;
    props = regionprops(imtranslate(mask,[0 rowSpacing/2]),'Centroid');
    props = props(~badSpots);
    for i=1:numel(props)
        bkgMask(round(props(i).Centroid(2))-2:round(props(i).Centroid(2))+2,round(props(i).Centroid(1))-2:round(props(i).Centroid(1))+2)=true;
    end
else
    bkgMask = [];
    computeBkg = false;
end
%% Preprocess Cy3 images and extract spot intensities
if computeBkg
    bkg = [];
end
[ValuesCy5, currBkg] = computeSpotIntensitiesFromRawData(dirCy5,parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,[],[],parameters.angularDisplacement,parameters.roi,offset,mask,bkgMask);
if computeBkg
   bkg = [bkg currBkg];
end
