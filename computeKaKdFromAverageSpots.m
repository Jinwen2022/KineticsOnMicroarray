function [ka,kd,eq,dna,dnaTable,kaStderror,kdStderror,eqStderror,countedspots]=computeKaKdFromAverageSpots(varargin)
% Computes average association and disassociation rates per operator
% sequence using the average spot method and linear fit.
%
% Input:
%   dnaSequences - 2D cell array with DNA sequences in microarray spots.
%   spotIntensities - 2D cell array, each element is a vector with average
%                     spot intensities of the corresponding microarray spot
%                     over acquisition frames.
%   operatorSequences - cell array with operator sequences
%   randSequence - char array with a random DNA sequences used for signal
%                  subtraction. If '', skip subtraction.
%   mutations - integer array with number of mutations
%   associationFramesLinearFit - frames for association rate linear fit
%   dissociationFramesLinearFit - frames for dissociation rate linear fit.
%                                 If empty, use only linear fit.
%   equilibirumFrames - frames for estimating spot intensity at equilibrium 
%   timeArray - vector with frame acqusition times
%   lacIConcentration - scalar
%   kdThreshold - discard sequences that are below kdThreshold at the
%                 initial dissociation frame. Optional, by default -Inf.
%
% Output:
%   ka - array with association rate values
%   kd - array with dissociation rate  values
%   eq - array with average operator intensity at equilibrium
%   dna - array with corresponding dna sequences
%   dnaTable - 3D logical array. dnaTable(i,j,k)=true if dna(i) differs
%              from operatorSequences(j) by mutations(k) gens.
%   kaStderror - standard error of mean of association rate estimates
%   kdStderror - standard error of mean of dissociation rate estimates
%   eqStderror - standard error of mean of equillibrium value estimates
%   countedspots - number of spots per analysed dna sequence.
%
% See also: computeKd, computeKa.
%% Parse parameters
ip = inputParser();
ip.addRequired('dnaSequences',@(x) iscell(x));
ip.addRequired('spotIntensities',@(x) isnumeric(x));
ip.addRequired('operatorSequences',@(x) iscell(x));
ip.addRequired('randSequence',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addRequired('mutations',@(x) isnumeric(x));
ip.addRequired('associationFramesLinearFit',@(x) isnumeric(x));
ip.addRequired('dissociationFramesLinearFit',@(x) isnumeric(x));
ip.addRequired('equilibriumFrames',@(x) isnumeric(x));
ip.addRequired('timeArray',@(x) isnumeric(x));
ip.addRequired('lacIConcentration',@(x) isnumeric(x) && numel(x)==1);
ip.addOptional('kdThreshold',-Inf,@(x) isnumeric(x) && numel(x)==1);
ip.addOptional('trimCutoff',30,@(x) isnumeric(x) && numel(x)==1);
ip.parse(varargin{:});

dnaSequences = ip.Results.dnaSequences;
spotIntensities = ip.Results.spotIntensities;
operatorSequences = ip.Results.operatorSequences;
randSeq = ip.Results.randSequence;
mutations = ip.Results.mutations;
framesALinearFit = ip.Results.associationFramesLinearFit;
framesDLinearFit = ip.Results.dissociationFramesLinearFit;
equilibriumFrames = ip.Results.equilibriumFrames;
time = ip.Results.timeArray;
lacIConcentration = ip.Results.lacIConcentration;
kdThreshold = ip.Results.kdThreshold;
trimCutoff = ip.Results.trimCutoff;
%% Compute average spot intensity per frame of the random sequence 
if isempty(randSeq)
    meanRandSeqValues = 0;
else
    randSeqIndex = strcmp(dnaSequences,randSeq);
    randSeqArrays = spotIntensities(:,randSeqIndex);
    meanRandSeqValues = trimmean(randSeqArrays,30,2);
%     figure, plot(randSeqArrays),hold on
%     plot(meanRandSeqValues,'k-','LineWidth',2)
end
%% Compute kd and ka*lacIConcentation using linear fit
maxMutationNum = length(operatorSequences{1});
[uDnaSequences,~,iu] = unique(dnaSequences);
dna = uDnaSequences;
dnaTable = false(numel(uDnaSequences),numel(operatorSequences),numel(mutations));
ka = nan(1,numel(uDnaSequences));
kd = nan(1,numel(uDnaSequences));
kaStderror = nan(1,numel(uDnaSequences));
kdStderror = nan(1,numel(uDnaSequences));
countedspots = nan(1,numel(uDnaSequences));
eq = nan(1,numel(uDnaSequences));
eqStderror = nan(1,numel(uDnaSequences));

for i=1:numel(uDnaSequences)
    if numel(uDnaSequences{i}) == maxMutationNum
        for j=1:numel(operatorSequences)
            currMutationNum = sum(uDnaSequences{i}~=operatorSequences{j});
            dnaTable(i,j,currMutationNum==mutations) = true;
        end
        if any(any(dnaTable(i,:,:)))
            seqIndex = iu==i;
            operatorValues = spotIntensities(:,seqIndex);
            ka_spots = nan(1,size(operatorValues,1));
            kd_spots = nan(1,size(operatorValues,1));
            eq_values = nan(1,size(operatorValues,1));
            for k = 1:size(operatorValues,2)
                %compute ka and kd using linear fit
                spotValues = operatorValues(:,k) - meanRandSeqValues;
                eq_values(k) = mean(spotValues(equilibriumFrames));
                ka_spots(k) = computeKa(spotValues,time,framesALinearFit,1);
                kd_spots(k) = computeKd(spotValues,time,framesDLinearFit,'linear','Threshold',kdThreshold);
            end
            index = eq_values>-Inf;
            ka_spots = ka_spots(index);
            kd_spots = kd_spots(index);
            eq_values = eq_values(index);
            [~,index1] = rmoutliers(ka_spots,'median');
            [~,index2] = rmoutliers(kd_spots,'median');
            [~,index3] = rmoutliers(eq_values,'median');
            ka_spots = ka_spots(~index1 & ~index2 & ~index3);
            kd_spots = kd_spots(~index1 & ~index2 & ~index3);
            eq_values = eq_values(~index1 & ~index2 & ~index3);
            ka(i) = mean(ka_spots);
            kd(i) = mean(kd_spots);
            eq(i) = mean(eq_values);
            kaStderror(i) = std(ka_spots)/sqrt(length(ka_spots));
            kdStderror(i) = std(kd_spots)/sqrt(length(kd_spots));
            eqStderror(i) = std(eq_values)/sqrt(length(eq_values));
            countedspots(i) = length(ka_spots);
        end
    end
end
ind = ~isnan(ka) & ~isnan(kd);
dna = dna(ind);
dnaTable = dnaTable(ind,:,:);
ka = ka(ind);
kd = kd(ind);
kaStderror = kaStderror(ind);
kdStderror = kdStderror(ind);
eq = eq(ind);
eqStderror = eqStderror(ind);
countedspots = countedspots(ind);
%% Get ka by dividing by lacI concentration
ka = ka/lacIConcentration;
kaStderror = kaStderror/lacIConcentration;