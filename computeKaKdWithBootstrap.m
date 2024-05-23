function [ka,kd,eq,dna,dnaTable,kaStderror,kdStderror,eqStderror,nSpots]=computeKaKdWithBootstrap(varargin)
% Computes average association and disassociation rates per operator
% sequence using the average curve method and linear fit. Standard errors
% of the mean are estimated using 500 bootstraping resamples.
%
% Input:
%   dnaSequences - 2D cell array with DNA sequences in microarray spots.
%   spotIntensities - 3D array, where (f,i,j) is the intensity of
%                     microarray spot (i,j) at acquisition frame f.
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
%   trimCutoff - compute trimmead average spot intensity per frame by
%                discarding trimcutoff/2 % highest and lowest data.
%                Optional, by default 30.
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
framesA = ip.Results.associationFramesLinearFit;
framesD = ip.Results.dissociationFramesLinearFit;
framesE = ip.Results.equilibriumFrames;
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
    meanRandSeqValues = trimmean(randSeqArrays,trimCutoff,2);
end
%% Compute kd and ka*lacIConcentation using linear fit
maxMutationNum = length(operatorSequences{1});
[uDnaSequences,~,iu] = unique(dnaSequences);
dna = uDnaSequences;
dnaTable = false(numel(uDnaSequences),numel(operatorSequences),numel(mutations));

%meanOperatorValues = cell(numel(uDnaSequences),1);
seqIndices = cell(numel(uDnaSequences),1);
for i=1:numel(uDnaSequences)
    if numel(uDnaSequences{i}) == maxMutationNum
        for j=1:numel(operatorSequences)
            currMutationNum = sum(uDnaSequences{i}~=operatorSequences{j});
            dnaTable(i,j,currMutationNum==mutations) = true;
        end
        if any(any(dnaTable(i,:,:)))
            seqIndices{i} = find(iu==i);
        end
    end
end
ind = ~cellfun('isempty',seqIndices);
dna = dna(ind);
dnaTable = dnaTable(ind,:,:);
seqIndices = seqIndices(ind);
%% Compute kd and ka*lacIConcentation using linear fit
ka = nan(1,numel(seqIndices));
kd = nan(1,numel(seqIndices));
eq = nan(1,numel(seqIndices));
kaStderror = nan(1,numel(seqIndices));
kdStderror = nan(1,numel(seqIndices));
eqStderror = nan(1,numel(seqIndices));
nSpots  = nan(1,numel(seqIndices));
rng(1);
Nboot = 500;
% create random number stream for reprodicibility
sc = parallel.pool.Constant(RandStream('Threefry'));
parfor i=1:numel(seqIndices)
    operatorValues = spotIntensities(:,seqIndices{i})-meanRandSeqValues;
    meanOperatorValues = trimmean(operatorValues,trimCutoff,2);
    ka(i) = computeKa(meanOperatorValues,time,framesA,1);
    kd(i) = computeKd(meanOperatorValues,time,framesD,'linear','Threshold',kdThreshold);
    eq(i) = mean(meanOperatorValues(framesE));
    %Bootstrap
    nSpots(i) = size(operatorValues,1);
    kaBootstr = nan(1,Nboot);
    kdBootstr = nan(1,Nboot);
    eqBootstr = nan(1,Nboot);
    stream = sc.Value;        % Extract the stream from the Constant
    stream.Substream = i;
    for j=1:Nboot
        spotIndices = ceil(rand(stream,1,nSpots(i))*nSpots(i));
        meanOperatorValues = trimmean(operatorValues(spotIndices,:),trimCutoff,2);
        kaBootstr(j) = computeKa(meanOperatorValues,time,framesA,1);
        kdBootstr(j) = computeKd(meanOperatorValues,time,framesD,'linear','Threshold',kdThreshold);
        eqBootstr(j) = mean(meanOperatorValues(framesE));
    end
    kaStderror(i) = std(rmoutliers(kaBootstr,'median'));
    kdStderror(i) = std(rmoutliers(kdBootstr,'median'));
    eqStderror(i) = std(rmoutliers(eqBootstr,'median'));
end
%% Remove nan values
ind = ~isnan(ka) & ~isnan(kd);
ka = ka(ind);
kd = kd(ind);
eq = eq(ind);
dna = dna(ind);
dnaTable = dnaTable(ind,:,:);
kaStderror = kaStderror(ind);
kdStderror = kdStderror(ind);
eqStderror = eqStderror(ind);
nSpots = nSpots(ind);
%% Get ka by deleting by lacI concentration
ka = ka/lacIConcentration;
kaStderror = kaStderror/lacIConcentration;