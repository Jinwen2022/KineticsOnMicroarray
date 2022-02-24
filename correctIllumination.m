function [weights, blurImage, averImage] = correctIllumination(location,gaussSigma,posRange,imRange)
% Computes pixel weights for correction of uneven illumination. For that,
% the function computes an average fluorescent image and smoothes it with a
% Gaussian filter. 
% The weights are computed as mean(blurImage)./blurImage
%
% Input:
%   location - path to folders with fluorescent images
%   gaussSigma - Gaussian smoothing kernel sigma
%   posRange - optional, range of positions (tiles). Default [], i.e. all
%              positions are used.
%   imRange - optional, range of images. Default [], i.e. all images are
%             used.
%
% Output:
%   weights - double-precision matrix with pixel weights, the same size as
%             microarray tiles.
%   blurImage - smoothened average image.
%   averImage - average image.

if nargin<3
    posRange = [];
end
if nargin<4
    imRange = [];
end

try
    imAdapterObj = genericReadAsFrames('metadata.txt',location);
    posList = imAdapterObj.getPositionList();
    if ~isempty(posRange)
        posList = posList(posRange);
    end
    chanNames = imAdapterObj.getChannels();
    chanName = chanNames{1};
    if isempty(imRange)
        nImagesPerPosition = zeros(size(posList));
        for i=1:length(posList)
            nImagesPerPosition(i) = numel(imAdapterObj.getIndicesByChanName(i,chanName));
        end
        imRange = 1:min(nImagesPerPosition);
    end
catch % metadata file is missing for some reasons..
    imds = imageDatastore(location,'IncludeSubfolders',true);
    [~,posList0] = fileparts(fileparts(fileparts(imds.Files)));
    posList = unique(posList0);
    [~, chanName] = fileparts(fileparts(imds.Files{1}));
    if ~isempty(posRange)
        posList = posList(posRange);
    end
    if isempty(imRange)
        nImagesPerPosition = zeros(size(posList));
        for i=1:length(posList)
            nImagesPerPosition(i) =  sum(strcmp(posList0,posList{i}));
        end
        imRange = 1:min(nImagesPerPosition);
    end
end

imLists = cell(size(posList));
for i=1:numel(posList)
    imDir = fullfile(location,posList{i},chanName);
    imList = dir(fullfile(imDir,'*.tif*'));
    imLists{i} = fullfile(imDir,{imList(imRange).name});
end
imLists = [imLists{:}];
imInfo = imfinfo(imLists{1});
w = imInfo(1).Height;
h = imInfo(1).Width;
ims = uint16(zeros(h,w,numel(posList)*numel(imRange)));
parfor i=1:numel(imLists)
    ims(:,:,i) = tiffread(imLists{i});
end

averImage = trimmean(ims,30,3);
blurImage = imgaussfilt(averImage,gaussSigma);
weights = mean(mean(blurImage))./blurImage;