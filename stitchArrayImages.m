function [mergedImages, stitchROI] = stitchArrayImages(location,overlap,rotationAngle,posColRange,imRange,posRot90,weights,offset)
% Stitches together MxN partially overlapping images.
% Images might have to be rotated a bit for properly alignment.
% Code parallelization is used when there are multiple frames.
%
% Input:
%   location - char array, path to microscopy data folder
%   overlap - even number, image overlap along x- and y-axis in pixels
%   rotationAngle - scalar, rotate images counterclockwise by
%                   rotationAngle degrees.
%   posColRange - column range of positions. Default [], use all positions
%   imRange - range of frames. Defaut [], use all frames
%   posRot90 - scalar, rotate 2D cell array with position list by
%              posRot90*90 degrees counterclockwise. Default -1.
%   weights - matrix with illumination correction pixel weights. Must be of
%             the same size as microscopy images.
%
% Output:
%   mergedImages - cell array with merged images if there are multiple
%                  frames, or a uint16 merged image.
%   stitchROI - MxNx8 array with border coordinates of MxN tiles:
%               [xGlobMin xGlobMax yGlobMin yGlobMax
%                xLocMin  xLocMax  yLocMin  yLocMax],
%               where "global" coordinates refer to the output images and
%               "local" coordinates refer to the raw images. Optional.

% Parse input
if nargin<4
    posColRange = [];
end
if nargin<5
    imRange = [];
end
if nargin<6
    posRot90 = -1;
end
if nargin<7
    weights = 1;
end
if nargin<8
    offset = 0;
end
if numel(weights)>1
    evenIllumination = true;
else
    evenIllumination = false;
end

if mod(overlap,2)==1
    error('Overlap must be an even number')
end
overlap2 = overlap/2;
try
%if exist(fullfile(location,'metadata.txt'),'file')
    imAdapterObj = genericReadAsFrames('metadata.txt',location);
    posList = imAdapterObj.getPositionList();
    chanNames = imAdapterObj.getChannels();
    chanName = chanNames{1};
    if isempty(imRange)
        nImagesPerPosition = zeros(size(posList));
        for i=1:length(posList)
%             imDir = fullfile(location,posList{i},chanName);
%             imList = dir(fullfile(imDir,'*.tif*'));
%             nImagesPerPosition(i) = numel(imList);
            nImagesPerPosition(i) = numel(imAdapterObj.getIndicesByPosAndChanName(posList{i},chanName));
        end
        imRange = 1:min(nImagesPerPosition);
    end
catch % metadata file is missing for some reasons..
    imds = imageDatastore(location,'IncludeSubfolders',true);
    [~,posList0] = fileparts(fileparts(fileparts(imds.Files)));
    posList = unique(posList0);
    [~, chanName] = fileparts(fileparts(imds.Files{1}));
    if isempty(imRange)
        nImagesPerPosition = zeros(size(posList));
        for i=1:length(posList)
            nImagesPerPosition(i) =  sum(strcmp(posList0,posList{i}));
        end
        imRange = 1:min(nImagesPerPosition);
    end
end
rows = str2double(cellfun(@(x) x(end-2:end),posList,'UniformOutput',false))+1;
cols = str2double(cellfun(@(x) x(end-6:end-4),posList,'UniformOutput',false))+1;
nRows = numel(unique(rows(~isnan(rows))));
nCols = numel(unique(cols(~isnan(cols))));
posList = reshape(posList,nRows,nCols);
posList = rot90(posList,posRot90);
if ~isempty(posColRange)
    posList = posList(:,posColRange);
end
[nRows, nCols] = size(posList);

imLists = cell(size(posList));
for i=1:numel(posList)
    imDir = fullfile(location,posList{i},chanName);
    imList = dir(fullfile(imDir,'*.tif*'));
    imLists{i} = fullfile(imDir,{imList(imRange).name});
end
imInfo = imfinfo(imLists{1}{1});
w = imInfo(1).Height;
h = imInfo(1).Width;
H = h*nRows -overlap*(nRows-1);
W = w*nCols -overlap*(nCols-1);

if numel(imRange)>1
    imListsPar = cell(1,numel(imRange));
    for k=1:numel(imRange)
        imListsPar{k} = cell(nRows,nCols);
        for i=1:nRows
            for j=1:nCols
                imListsPar{k}(i,j) = imLists{i,j}(k);
            end
        end
    end
    mergedImages = cell(1,numel(imRange));
    parfor k=1:numel(imRange)
        currMergedImage = uint16(zeros(H,W));
        for i=1:nRows
            if i==1
                rowInd = 1:h-overlap2;
                subImRowInd = 1:h-overlap2;
            elseif i==nRows
                rowInd = (i-1)*(h-overlap)+overlap2+1:H;
                subImRowInd = overlap2+1:h;
            else
                rowInd = (i-1)*(h-overlap)+overlap2+1:i*(h-overlap)+overlap2;
                subImRowInd = overlap2+1:h-overlap2;
            end
            for j=1:nCols
                if j==1
                    colInd = 1:w-5;%-overlap2;
                    subImColInd = 1:w-5;%-overlap2;
                elseif j==nCols
                    colInd = (j-1)*(w-overlap)+overlap-5+1:W;
                    subImColInd = overlap-5+1:w;
                else
                    colInd =(j-1)*(w-overlap)+overlap-5+1:j*(w-overlap)+overlap-5;
                    subImColInd = overlap-5+1:w-5;%-overlap2;
                end
                im = tiffread(imListsPar{k}{i,j});
                im = im-offset;
                if rotationAngle
                    im = imrotate(im,rotationAngle,'bilinear','crop');
                end
                if evenIllumination
                    im = uint16(double(im).*weights);
                end
                currMergedImage(rowInd,colInd) = im(subImRowInd,subImColInd);
            end
        end
        mergedImages{k} = currMergedImage;
    end
else
    mergedImages = uint16(zeros(H,W));
    for i=1:nRows
        if i==1
            rowInd = 1:h-overlap2;
            subImRowInd = 1:h-overlap2;
        elseif i==nRows
            rowInd = (i-1)*(h-overlap)+overlap2+1:H;
            subImRowInd = overlap2+1:h;
        else
            rowInd = (i-1)*(h-overlap)+overlap2+1:i*(h-overlap)+overlap2;
            subImRowInd = overlap2+1:h-overlap2;
        end
        for j=1:nCols
            if j==1
                colInd = 1:w-5;%-overlap2;
                subImColInd = 1:w-5;%-overlap2;
            elseif j==nCols
                colInd = (j-1)*(w-overlap)+overlap-5+1:W;
                subImColInd = overlap-5+1:w;
            else
                colInd =(j-1)*(w-overlap)+overlap-5+1:j*(w-overlap)+overlap-5;
                subImColInd = overlap-5+1:w-5;%-overlap2;
            end
            im = tiffread(imLists{i,j}{1});
            im = im-offset;
            if rotationAngle
                im = imrotate(im,rotationAngle,'bilinear','crop');
            end
            if evenIllumination
                im = uint16(double(im).*weights);
            end
            mergedImages(rowInd,colInd) = im(subImRowInd,subImColInd);
        end
    end
end
if nargout==2
     stitchROI = zeros(nRows,nCols,8);
     for i=1:nRows
        if i==1
            rowInd = 1:h-overlap2;
            subImRowInd = 1:h-overlap2;
        elseif i==nRows
            rowInd = (i-1)*(h-overlap)+overlap2+1:H;
            subImRowInd = overlap2+1:h;
        else
            rowInd = (i-1)*(h-overlap)+overlap2+1:i*(h-overlap)+overlap2;
            subImRowInd = overlap2+1:h-overlap2;
        end
        for j=1:nCols
            if j==1
                colInd = 1:w-5;%-overlap2;
                subImColInd = 1:w-5;%-overlap2;
            elseif j==nCols
                colInd = (j-1)*(w-overlap)+overlap-5+1:W;
                subImColInd = overlap-5+1:w;
            else
                colInd =(j-1)*(w-overlap)+overlap-5+1:j*(w-overlap)+overlap-5;
                subImColInd = overlap-5+1:w-5;%-overlap2;
            end
            stitchROI(i,j,:)=[colInd(1) colInd(end) rowInd(1) rowInd(end)...
                              subImColInd(1) subImColInd(end) subImRowInd(1) subImRowInd(end)];
        end
     end
end
