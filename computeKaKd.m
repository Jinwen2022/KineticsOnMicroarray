function [ka,kd,dna,dnaTable,kaInit,kdInit]=computeKaKd(varargin)
% Computes average association and disassociation rates per operator
% sequence. For each operator, it computes an average spot intesity curve
% in each frame. Initial ka and kd values are found using linear fit, they
% optionally can be further refined using exponential fit.
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
%   associationFramesExpFit - frames for association rate exponential fit
%   dissociationFramesLinearFit - frames for dissociation rate linear fit.
%                                 If empty, use only linear fit.
%   dissociationFramesExpFit - frames for dissociation rate exponential
%                              fit. If empty, use only linear fit.
%   timeArray - vector with frame acqusition times
%   lacIConcentration - scalar
%   kdThreshold - discard sequences that are below kdThreshold at the
%                 initial dissociation frame.
%
% Output:
%   ka - array with association rate values
%   kd - array with dissociation rate  values
%   dna - array with corresponding dna sequences
%   dnaTable - 3D logical array. dnaTable(i,j,k)=true if dna(i) differs
%              from operatorSequences(j) by mutations(k) gens.
%   kaInit - array with association rate values using linear fit
%   kdInit - array with dissociation rate  values using linear fit.
%
% See also: computeKd, associationCurve
%% Parse parameters
ip = inputParser();
ip.addRequired('dnaSequences',@(x) iscell(x));
ip.addRequired('spotIntensities',@(x) iscell(x));
ip.addRequired('operatorSequences',@(x) iscell(x));
ip.addRequired('randSequence',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addRequired('mutations',@(x) isnumeric(x));
ip.addRequired('associationFramesLinearFit',@(x) isnumeric(x));
ip.addRequired('associationFramesExpFit',@(x) isnumeric(x));
ip.addRequired('dissociationFramesLinearFit',@(x) isnumeric(x));
ip.addRequired('dissociationFramesExpFit',@(x) isnumeric(x));
ip.addRequired('timeArray',@(x) isnumeric(x));
ip.addRequired('lacIConcentration',@(x) isnumeric(x) && numel(x)==1);
ip.addRequired('kdThreshold',@(x) isnumeric(x) && numel(x)==1);
ip.parse(varargin{:});

dnaSequences = ip.Results.dnaSequences;
spotIntensities = ip.Results.spotIntensities;
operatorSequences = ip.Results.operatorSequences;
randSeq = ip.Results.randSequence;
mutations = ip.Results.mutations;
framesALinearFit = ip.Results.associationFramesLinearFit;
framesAExpFit = ip.Results.associationFramesExpFit;
framesDLinearFit = ip.Results.dissociationFramesLinearFit;
framesDExpFit = ip.Results.dissociationFramesExpFit;
time = ip.Results.timeArray;
lacIConcentration = ip.Results.lacIConcentration;
kdThreshold = ip.Results.kdThreshold;
%% Compute average spot intensity per frame of the random sequence 
if isempty(randSeq)
    meanRandSeqValues = 0;
else
    randSeqIndex = strcmp(dnaSequences,randSeq);
    randSeqArrays = cell2mat(spotIntensities(randSeqIndex));
    meanRandSeqValues = trimmean(randSeqArrays,30,1);
end
%% Compute kd and ka*lacIConcentation using linear fit
maxMutationNum = length(operatorSequences{1});
[uDnaSequences,~,iu] = unique(dnaSequences);
dna = uDnaSequences;
dnaTable = false(numel(uDnaSequences),numel(operatorSequences),numel(mutations));

meanOperatorValues = cell(numel(uDnaSequences),1);
for i=1:numel(uDnaSequences)
    if numel(uDnaSequences{i}) == maxMutationNum
        for j=1:numel(operatorSequences)
            currMutationNum = sum(uDnaSequences{i}~=operatorSequences{j});
            dnaTable(i,j,currMutationNum==mutations) = true;
        end
        if any(any(dnaTable(i,:,:)))
            seqIndex = iu==i;
            operatorValues = cell2mat(spotIntensities(seqIndex));
            meanOperatorValues{i} = trimmean(operatorValues,30,1);
            meanOperatorValues{i} = meanOperatorValues{i}-meanRandSeqValues;
        end
    end
end
ind = ~cellfun('isempty',meanOperatorValues);
dna = dna(ind);
dnaTable = dnaTable(ind,:,:);
meanOperatorValues = meanOperatorValues(ind);
if ~isempty(framesAExpFit)
    fArray = cellfun(@(x) x(framesAExpFit),meanOperatorValues,'UniformOutput',false);
end
%% Compute kd and ka*lacIConcentation using linear fit
ka = nan(1,numel(uDnaSequences));
kd = nan(1,numel(uDnaSequences));
parfor i=1:numel(meanOperatorValues)
    ka(i) = computeKa(meanOperatorValues{i},time,framesALinearFit,1);
    kd(i) = computeKd(meanOperatorValues{i},time,framesDLinearFit,'linear','Threshold',kdThreshold);
end
%% Compute kd using exp fit
kdInit = kd;
if ~isempty(framesDExpFit)
    parfor i=1:numel(meanOperatorValues)
        kd(i) = computeKd(meanOperatorValues{i},time,framesDExpFit,'exp','Threshold',kdThreshold,'StartPoint',kd(i));
    end
end
%% Remove nan values
ind = ~isnan(ka) & ~isnan(kd);
ka = ka(ind);
kd = kd(ind);
kdInit = kdInit(ind);
dna = dna(ind);
dnaTable = dnaTable(ind,:,:);
meanOperatorValues = meanOperatorValues(ind);
%% Compute ka using exp fit
kaInit = ka;
if ~isempty(framesAExpFit)
    fArray = fArray(ind);
    F = cell2mat(fArray');
    t = time(framesAExpFit)-time(framesAExpFit(1));
    nt = numel(t);
    nka = numel(ka);
    % Set optimization parameters
    options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
    tolF = 1e-3; 
    tolX = 1e-6;
    maxIter = 10000;

    % Compute a init
    Ffit = nan(size(F));
    a = 1;
    for i=1:nka
        ind = (i-1)*nt+1:i*nt;
        Ffit(ind) = associationCurve(ka(i),a,kd(i),t,0);
    end
    a0 = a;
    a = F/Ffit*a0;
    Ffit = Ffit*(a/a0);
    res=norm(F-Ffit);

    iter=0;
    dx=inf;
    fprintf('iter\tres\t\tstep\n');
    fprintf('%d\t%.4e\t\n',iter,res);
    tic;
    while res>tolF && iter<=maxIter && dx>tolX
        iter= iter+1;
        kaArray0 = ka;
        parfor i=1:nka
            ind = (i-1)*nt+1:i*nt;
            fun = @(x) associationCurve(x,a,kd(i),t,F(ind));
            ka0 = ka(i);
            ka(i) = lsqnonlin(fun,ka0,0,[],options);
        end
        for i=1:nka
            ind = (i-1)*nt+1:i*nt;
            Ffit(ind) = associationCurve(ka(i),a,kd(i),t,0);
        end
        a0 = a;
        a = F/Ffit*a0;
        Ffit = Ffit*(a/a0);
        res = norm(F-Ffit);
        dx = norm(ka-kaArray0);
        fprintf('%d\t%.4e\t%.4e\n',iter,res,dx);
    end
    toc;
end
%% Get ka by deleting by lacI concentration
ka = ka/lacIConcentration;
kaInit = kaInit/lacIConcentration;
%% Plot fitted curves for the selected operators with 0 mutations
col = find(mutations==0, 1);
if ~isempty(col) && ~isempty(framesAExpFit)
    figure
    tiledlayout(numel(operatorSequences),2);
    for k=1:numel(operatorSequences)
        i = find(strcmp(dna,operatorSequences{k}));
        ax=nexttile(2*k-1);
        if ~isempty(framesDExpFit)
            computeKd(meanOperatorValues{i},time,framesDExpFit,'exp','Threshold',kdThreshold,'StartPoint',kdInit(i),'Axes',ax,'MakePlot',true);
        else
            computeKd(meanOperatorValues{i},time,framesDLinearFit,'linear','Threshold',kdThreshold,'Axes',ax,'MakePlot',true);
        end
        ax = nexttile(2*k);
        hold(ax,'on');
        ind = (i-1)*nt+1:i*nt;
        plot(t,F(ind),'o',t,Ffit(ind))
%         [~, p] = computeKa(meanOperatorValues{i},time,framesALinearFit,1,false);
        F0 = associationCurve(kaInit(i)*lacIConcentration,1,0,t,0);%polyval(p,time(framesAExpFit));
        plot(t,F0,'r--');
        hold(ax,'off');
    end
end