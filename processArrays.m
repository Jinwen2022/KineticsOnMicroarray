function [valuesCy5,valuesCy3,mask,dnaSequences,nLabeledNucleotides] = processArrays(dirCy5,dirCy3,dirOffset,dnaFile,outputDir,parameters)
% Performs analysis of specific operator binding experiments using
% microarrays. The function stores in the output analysis folder merged
% Cy3 and Cy5 images of the microarray, and file 'results.mat' with output
% variables.
%
% Input:
%   dirCy5 - path to directory with Cy5 fluorescence images. Leave empty if
%            no images is available.
%   dirCy3 - cell array, each element is a path to directory with Cy3
%            fluorescence images.
%   dirOffset - path to directory with an offset Cy3 image, used for
%               background subtraction. If empty no background subtraction.
%   dnaFile - text file with a table of dna sequences in the microarray, or
%             cell array with two files' paths: the first geo file is with
%             DNA keys and spot locations, and the second file with DNA
%             keys and DNA sequences.
%   outputDir - path to to directory where the results will be saved.
%   parameters - scructure array with fields:
%       overlap - number of tile overlapping pixels. If the overlap in Cy3
%                 and Cy5 is different, use overlapCy5 and overlapCy3.
%       overlapCy5 - number of tile overlapping pixels in Cy5. Optional.
%       overlapCy3 - number of tile overlapping pixels in Cy3. Optional.
%       rangePositionsCy5 - array with position indices in Cy5. If empty,
%                           all positions are used. Default [].
%       weightsCy5 - matrix for illumination correction in Cy5. If empty,
%                    skip correction. Default [].
%       angularDisplacementCy5Tile - in degrees, rotate counter-clockwise
%                                    raw Cy5 images (tiles).
%       angularDisplacementCy5 - in degrees, rotate counter-clockwise
%                                merged Cy5 image.
%       roiCy5 - cropping region in Cy5 image after rotation.
%       rangePositionsCy3 - array with position indices in Cy5. If empty,
%                           all positions are used. Default [].
%       weightsCy3 - matrix for illumination correction in Cy3.  If empty,
%                    skip correction. Default [].
%       angularDisplacementCy3Tile - in degrees, rotate counter-clockwise
%                                    raw Cy3 images (tiles).
%       angularDisplacementCy3 - in degrees, rotate counter-clockwise
%                                merged Cy3 images.
%       roiCy5 - cropping region in Cy3 images after rotation.
%       pxSize - pixel size in mum.
%       topLeftCorner - [x y] coordinates of the centroid of the
%                       (imaginary) spot in the top left corner of
%                       preprocessed image.
%       labeledNucleotide - char, labeled nucleotide.
%       dnaLength - length of valid dna sequences. Default 60.
%       primer - char array, remove the primer from dna sequences to count
%                the labeled nucleotide. Default 'GTCTGTGTTCCGTTGTCCGTGCTG'
%
% Output:
%   valuesCy5 - matrix with microarray average spot intensities in Cy5.
%   valuesCy3 - cell 2D array, each element corresponds to a microarray
%               spot and contains a vector of spot intensities in Cy3
%               sorted by image acquisition time.
%   mask - boolean matrix of microarray spots in the preprocessed images.
%   dnaSequences - cell array with dna sequences in the microarray spots.
%   nLabeledNucleotides - number of labeled nucleotide counts per spot.
%
% Dependencies: calls some functions from ImAnalysis.
%
% Spartak Zikrin, Elf lab, 2021-06-02.
%% Parse parameters
if ~isfield(parameters,'overlapCy5')
    parameters.overlapCy5 = parameters.overlap;
end
if ~isfield(parameters,'overlapCy3')
    parameters.overlapCy3 = parameters.overlap;
end
if ~isfield(parameters,'rangePositionsCy5')
    parameters.rangePositionsCy5 = [];
end
if ~isfield(parameters,'rangePositionsCy3')
    parameters.rangePositionsCy3 = [];
end
if ~isfield(parameters,'weightsCy5')
    parameters.weightsCy5 = [];
end
if ~isfield(parameters,'weightsCy3')
    parameters.weightsCy3 = [];
end
%% Preprocess Cy5 images
if ~isempty(dirCy5)
    t0=tic;
    slideCy5 = stitchArrayImages(dirCy5,parameters.overlapCy5,parameters.angularDisplacementCy5Tile,parameters.rangePositionsCy5,1,-1,parameters.weightsCy5);
    slideCy5 = preprocessImage(slideCy5,parameters.angularDisplacementCy5,parameters.roiCy5);
    imwrite(slideCy5,fullfile(outputDir,'cy5.tiff'));
    % Get the spot grid and compute Cy5 spot intensities
    mask = arrayGridMask(slideCy5,parameters.pxSize,parameters.topLeftCorner);
    valuesCy5 = computeSpotIntensities(slideCy5,mask);
    t=toc(t0);
    fprintf('Processed Cy5 image in %.3f s\n',t);
else
    slideCy5 = zeros(parameters.roiCy3([4 3])+1);
    mask = arrayGridMask(slideCy5,parameters.pxSize,parameters.topLeftCorner);
    valuesCy5 = [];
end
%% Preprocess offset image if it is available
t0=tic;
if ~isempty(dirOffset)
    offset = stitchArrayImages(dirOffset,parameters.overlapCy3,parameters.angularDisplacementCy3Tile,parameters.rangePositionsCy3,1,-1,parameters.weightsCy3);
    offset = preprocessImage(offset,parameters.angularDisplacementCy3,parameters.roiCy3);
    imwrite(offset,fullfile(outputDir,'offset.tiff'));
else
    offset = 0;
end
t=toc(t0);
fprintf('Processed offset image in %.3f s\n',t);
%% Preprocess Cy3 images together
t0=tic;
mkdir(fullfile(outputDir,'Cy3'));
angularDisplacementCy3 = parameters.angularDisplacementCy3;
roiCy3 = parameters.roiCy3;
for i=1:numel(dirCy3)
    currSlidesCy3 = stitchArrayImages(dirCy3{i},parameters.overlapCy3,parameters.angularDisplacementCy3Tile,parameters.rangePositionsCy3,[],-1,parameters.weightsCy3);
    if ~iscell(currSlidesCy3)
        currSlidesCy3 = {currSlidesCy3};
    end
    parfor j=1:numel(currSlidesCy3)
        imFilename = fullfile(outputDir,'Cy3',sprintf('cy3_%03d_%03d.tiff',i,j));
        if ~exist(imFilename,'file') % avoid recomputing images
            slideCy3 = preprocessImage(currSlidesCy3{j},angularDisplacementCy3,roiCy3);
            slideCy3 = slideCy3-offset;
            imwrite(slideCy3,imFilename);
        end
    end
end
% Concatenate slides
t=toc(t0);
fprintf('Merged Cy3 images in %.3f s\n',t);
%% Get Cy3 values
t0=tic;
imDir = fullfile(outputDir,'Cy3');
imList = dir(fullfile(imDir,'*.tif*'));
imList = {imList.name};
tmpValuesCy3 = cell(1,length(imList));
parfor k = 1:length(imList)
    slideCy3 = tiffread(fullfile(imDir,imList{k}));
    tmpValuesCy3{k} = computeSpotIntensities(slideCy3,mask);
end

[nRows,nCols] = size(tmpValuesCy3{1});
valuesCy3 = cell(nRows, nCols);
for i = 1:nRows
    for j = 1:nCols
        valuesCy3{i,j} = zeros(1,length(tmpValuesCy3));
        for k = 1:length(tmpValuesCy3)
            valuesCy3{i,j}(k)= tmpValuesCy3{k}(i,j);
        end 
    end
end
t=toc(t0);
fprintf('Computed Cy3 spot intensities in %.3f s\n',t);
%% Import DNA sequences
status = true;
if (ischar(dnaFile) || isstring(dnaFile)) && exist(dnaFile,'file')
    dnaSequences = getDnaSequences(dnaFile);
elseif iscell(dnaFile) && numel(dnaFile)==2
    dnaSequences = getDnaSequences(dnaFile{1},dnaFile{2});
else
    disp('Something is wrong with provided DNA file(s)');
    status = false;
end
if status
    nLabeledNucleotides = getNumLabeledNucleotides(dnaSequences,parameters.labeledNucleotide,parameters.dnaLength,parameters.primer);
else
    nLabeledNucleotides = [];
end
%% Save results
save(fullfile(outputDir,'results.mat'),'valuesCy5','valuesCy3','mask','dnaSequences','nLabeledNucleotides','parameters')