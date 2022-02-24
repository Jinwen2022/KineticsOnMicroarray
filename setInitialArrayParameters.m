function [parameters, slideCy5, slideCy3, mask] = setInitialArrayParameters(dirCy5,dirCy3,outputDir,parameters)
% This is a wrapper to set microarray image preprocessing parameters 
% such as cropping region, rotation angle, and top left spot coordinates.
%
% Input:
%   dirCy5 - char array, path to Cy5 images
%   dirCy3 - cell array, paths to Cy3 images
%   outputDir - output directory path
%   parameters - structure array with pipeline parameters
%
% Output:
%   parameters - structure array with updated pipeline parameters
%   slideCy5 - merged Cy5 image
%   slideCy3 - merged Cy3 image
%   mask - binary mask of microarray spots

if ~isempty(dirCy5)
    slideCy5 = stitchArrayImages(dirCy5,parameters.overlap,parameters.angularDisplacementCy5Tile,[],1,parameters.posRot90);
    imwrite(slideCy5,fullfile(outputDir,'cy5.tiff'));
    [slideCy5, parameters.angularDisplacementCy5,parameters.roiCy5] = preprocessImage(slideCy5,parameters.angularDisplacementCy5,parameters.roiCy5);
else
    slideCy5 = [];
end

if ~isempty(dirCy3)
    fullSlideCy3 = stitchArrayImages(dirCy3,parameters.overlap,parameters.angularDisplacementCy3Tile,[],1,parameters.posRot90);
    [slideCy3,parameters.angularDisplacementCy3,parameters.roiCy3] = preprocessImage(fullSlideCy3,parameters.angularDisplacementCy3,parameters.roiCy3);
    if ~isempty(dirCy5)
        parameters.roiCy3(3:4) = parameters.roiCy5(3:4);
        slideCy3 = preprocessImage(fullSlideCy3,parameters.angularDisplacementCy3,parameters.roiCy3);
        imwrite(slideCy3,fullfile(outputDir,'cy3.tiff'));
    end
end

if ~isempty(dirCy5)
    [mask, parameters.topLeftCorner] = arrayGridMask(slideCy5,parameters.pxSize,parameters.topLeftCorner,1);
else
    [mask, parameters.topLeftCorner] = arrayGridMask(slideCy3,parameters.pxSize,parameters.topLeftCorner,1);
end
figure, imshow(imadjust(slideCy3)); hold on
B = bwboundaries(mask,'noholes');
visboundaries(B,'EnhanceVisibility',0);