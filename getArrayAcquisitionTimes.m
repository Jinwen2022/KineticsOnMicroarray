function arrayAcquisitionTimes = getArrayAcquisitionTimes(locations,posRange)
% Extract microarray image acquisition times. 
% If there are multiple positions (tiles), use images from the first one.
%
% Input:
%   locations - cell array of char arrays, or char array with paths to
%               image folders
%   posRange - optional, array with position indices. Default [].
%
% Output:
%   arrayAcquisitonTimes - array with image acquisition times in seconds.

if nargin<2
    posRange = [];
end
if ~iscell(locations)
    locations = {locations};
end

arrayAcquisitionTimes = cell(numel(locations),1);
for dirIdx=1:numel(locations)
    imAdapterObj = genericReadAsFrames('metadata.txt', locations{dirIdx});
    posList = imAdapterObj.getPositionList();
    chanNames = imAdapterObj.getChannels();
    chanName = chanNames{1};
    if ~isempty(posRange)
        posList = posList(posRange);
    end
    nImagesPerPosition = zeros(size(posList));
    for i=1:length(posList)
%         imDir = fullfile(locations{dirIdx},posList{i},chanName);
%         imList = dir(fullfile(imDir,'*.tif*'));
%         nImagesPerPosition(i) = numel(imList);
        nImagesPerPosition(i) = numel(imAdapterObj.getIndicesByChanName(i,chanName));
    end
    imRange = 1:min(nImagesPerPosition);
    arrayAcquisitionTimes{dirIdx} = imAdapterObj.getAcquireTimes(posList{1},chanName);
    arrayAcquisitionTimes{dirIdx} = arrayAcquisitionTimes{dirIdx}(imRange);
end
arrayAcquisitionTimes = cell2mat(arrayAcquisitionTimes)';
arrayAcquisitionTimes = (arrayAcquisitionTimes-arrayAcquisitionTimes(1))/1000;