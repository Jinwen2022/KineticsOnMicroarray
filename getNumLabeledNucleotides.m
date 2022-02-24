function nLabeledNucleotides = getNumLabeledNucleotides(dnaSequences,labeledNucleotide,universalDnaLength,cutPrimer)
% Count number of a labeled nucleotide appearances in the DNA sequences
%
% Input:
%   dnaSequences - cell 2D array with DNA sequences in microarray spots
%   labeledNucleotide - char symbol, nucleotide
%   universalDnaLength - scalar, defines number of char symbols in valid
%                        DNA sequences. Default 60.
%   cutPrimer - char array with a primer that is removed from the DNA
%               sequences before the counting.
%               Default 'GTCTGTGTTCCGTTGTCCGTGCTG'.

% Parse parameters and set defaults
if nargin<3 || isempty(universalDnaLength)
    universalDnaLength = 60;
end
if nargin<4
    cutPrimer = 'GTCTGTGTTCCGTTGTCCGTGCTG';
end

% Find valid DNA sequences
dnaLengths = cellfun('length',dnaSequences);
ind = dnaLengths == universalDnaLength;
validDnaSequences = dnaSequences(ind);
if ~isempty(cutPrimer)
    validDnaSequences = erase(validDnaSequences,cutPrimer);
end

% Count number of labeled nucleotide in the valid DNA sequences
nLabeledNucleotides = nan(size(dnaSequences));
nLabeledNucleotides(ind) = count(validDnaSequences,labeledNucleotide);