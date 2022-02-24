function dnaSequencesCameraView = getDnaSequences(varargin)
% Reads a table with DNA sequences on microarray and
% converts it into a 2D cell array with DNA sequences on 
% preprocessed microarray images.
%
% getDnaSequences(filename) parses *filename* file which contains a table
% of size either 96x164 or 266x170 with DNA sequences on the microarray.
%
% getDnaSequences(geoFilename,sequenceFilename) parses info from two files:
% geoFilename file has a table with columns [id row col dna_key flag];
% sequenceFilename file has a table with columns [dna_key dna_sequence].
% The function uses dna_key values to set dna_sequence values in the
% microarray spots defined by the corresponding row and col values.
%
% Output:
%   dnaSequencesCameraView - 2D cell array with DNA sequences

if numel(varargin)==1
    dnaSequencesSlideView = readcell(varargin{1});
    dnaSequencesCameraView = rot90(dnaSequencesSlideView,2)';
else
    geoList = table2cell(readtable(varargin{1},'Delimiter','tab'));
    geoList = geoList(~isnan(cell2mat(geoList(:,2))),:);
    dnaSequenceList = table2cell(readtable(varargin{2},'Delimiter','tab'));
    [~,IG,ID] = intersect(geoList(:,4),dnaSequenceList(:,1),'stable');
    if size(geoList,1) == 96*164
        dnaSequencesSlideView = cell(96,164);
    elseif size(geoList,1) == 266*170
        dnaSequencesSlideView = cell(266,170);
    else
        error('There must be 96x164 or 266x170 microarray spots')
    end
    IG2 = false(size(geoList,1),1);
    IG2(IG) = true;
    ind = sub2ind(size(dnaSequencesSlideView),cell2mat(geoList(IG2,2)),cell2mat(geoList(IG2,3)));
    dnaSequencesSlideView(ind) = dnaSequenceList(ID,2);
    ind = sub2ind(size(dnaSequencesSlideView),cell2mat(geoList(~IG2,2)),cell2mat(geoList(~IG2,3)));
    dnaSequencesSlideView(ind) = geoList(~IG2,4);
    dnaSequencesCameraView = rot90(dnaSequencesSlideView,2);
end