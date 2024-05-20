function valuesCy3 = computeSpotIntensitiesFromRawData(location,overlap,rotAngleTile,posColRange,imRange,weights,rotAngle,roi,offsetIm,mask,bkgMask)
% This function runs a parfor loop over fluorescent frames of a microarray
% in a directory. At each iteration, it merges tiles into a full image in
% the corresponding frame, preprocess the image and extracts average
% fluorescence intensity per microarray spot.
%
% Input:
%   location - char array, path to microscopy data folder
%   overlap - even number, image overlap along x- and y-axis in pixels
%   rotAngleTile - scalar, rotate images counterclockwise by rotAngleTile degrees.
%   posColRange - column range of positions. If [] use all positions.
%   imRange - range of frames. If [] use all frames.
%   weights - matrix with illumination correction pixel weights. 
%             If [], skip illumination correction.
%   rotAngle - scalar, defines counter-clockwise image rotation in degrees.
%   roi - 4-element vector [XMIN YMIN W H] for image cropping region.
%   offsetIm -uint16 array (H+1)x(W+1) for background
%             subrtraction. Or simply 0 scalar. 
%   mask - binary array (H+1)x(W+1) of microarray spots.
%   bkgMask - binary array (H+1)x(W+1) of background spots.
% 
% Output:
%   valuesCy3 - 3D array,  where (f,i,j) is the intensity of microarray
%               spot (i,j) at acquisition frame f.
%
% See also: stitchArrayImages, computeSpotIntensities, preprocessImage
%% Parse parameters

if mod(overlap,2)==1
    error('Overlap must be an even number')
end
overlap2 = overlap/2;

if isempty(weights)
    evenIllumination = false;
else
    evenIllumination = true;
end
%% Get image info
try
    imAdapterObj = genericReadAsFrames('metadata.txt',location);
    posList = imAdapterObj.getPositionList();
    chanNames = imAdapterObj.getChannels();
    chanName = chanNames{1};
    if isempty(imRange)
        nImagesPerPosition = zeros(size(posList));
        for i=1:length(posList)
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
posList = rot90(posList,-1);
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

imListsPar = cell(1,numel(imRange));
for k=1:numel(imRange)
    imListsPar{k} = cell(nRows,nCols);
    for i=1:nRows
        for j=1:nCols
            imListsPar{k}(i,j) = imLists{i,j}(k);
        end
    end
end
%% Main code
tmpValuesCy3 =  cell(1,numel(imRange));
parfor k=1:numel(imRange)
    mergedIm = uint16(zeros(H,W));
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
                colInd = 1:w-overlap2;
                subImColInd = 1:w-overlap2;
            elseif j==nCols
                colInd = (j-1)*(w-overlap)+overlap2+1:W;
                subImColInd = overlap2+1:w;
            else
                colInd =(j-1)*(w-overlap)+overlap2+1:j*(w-overlap)+overlap2;
                subImColInd = overlap2+1:w-overlap2;
            end
            im = tiffread(imListsPar{k}{i,j});
            if evenIllumination
                im = uint16(double(im).*weights);
            end
            if rotAngleTile
                im = imrotate(im,rotAngleTile,'bilinear','crop');
            end
            mergedIm(rowInd,colInd) = im(subImRowInd,subImColInd);
        end
    end
    mergedIm = preprocessImage(mergedIm,rotAngle,roi);
    mergedIm = mergedIm-offsetIm;
    tmpValuesCy3{k} = computeSpotIntensities(mergedIm,mask,bkgMask);
end
%% Rearrange spot intensities array
[nRows,nCols] = size(tmpValuesCy3{1});
valuesCy3 = zeros(numel(tmpValuesCy3),nRows, nCols);
for k = 1:length(tmpValuesCy3)
    for i = 1:nRows
        for j = 1:nCols
            valuesCy3(k,i,j)= tmpValuesCy3{k}(i,j);
        end 
    end
end