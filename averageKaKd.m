function [kaArray,kdArray,dnaArray,kaCVArray,kdCVArray]=averageKaKd(varargin)
%% Parse parameters
if length(varargin)>1 && isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = varargin{1};
    makePlot = 1;
    varargin = varargin(2:end);
else
    makePlot = 0;
end
ip = inputParser();
ip.addRequired('dnaSequences',@(x) iscell(x));
ip.addRequired('spotIntensities',@(x) iscell(x));
ip.addRequired('operatorSequence',@(x) ischar(x) || isstring(x));
ip.addRequired('randSequence',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addRequired('Osym_DNA',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addRequired('mutationNum',@(x) isnumeric(x) && numel(x)==1);
ip.addRequired('associationFrames',@(x) isnumeric(x));
ip.addRequired('disassociationFrames',@(x) isnumeric(x));
ip.addRequired('timeArray',@(x) isnumeric(x));
ip.addRequired('lacIConcentration',@(x) isnumeric(x) && numel(x)==1);
ip.addParameter('Color','', @(x) iscell(x) || ischar(x) || isnumeric(x));
ip.addParameter('Marker','o',@(x) ischar(x))
ip.addParameter('MakePlot',makePlot)
ip.parse(varargin{:});

dnaSequences = ip.Results.dnaSequences;
spotIntensities = ip.Results.spotIntensities;
operatorSequence = ip.Results.operatorSequence;
randSeq = ip.Results.randSequence;
Osym_DNA = ip.Results.Osym_DNA;
mutationNum = ip.Results.mutationNum;
associationFrames = ip.Results.associationFrames;
disassociationFrames = ip.Results.disassociationFrames;
timeArray = ip.Results.timeArray;
lacIConcentration = ip.Results.lacIConcentration;
color = ip.Results.Color;
marker = ip.Results.Marker;
makePlot = ip.Results.MakePlot;

if makePlot && ~exist('ax','var')
    figure, ax=gca;
end

%% Calculate ka and kd

if isempty(randSeq)
    meanRandSeqValue = 0;
else
    randSeqIndex = strcmp(dnaSequences,randSeq);
    randSeqArrays = cell2mat(spotIntensities(randSeqIndex));
    meanRandSeqValue = trimmean(randSeqArrays,30,1);
end

[uDnaSequences,~,iu] = unique(dnaSequences);


Osym_index = find(strcmp(uDnaSequences,Osym_DNA));
seqIndex = iu==Osym_index;
operatorValues = cell2mat(spotIntensities(seqIndex));
ka_Osym=nan(1,length(operatorValues(:,1)));
for j = 1:size(operatorValues,1)
    meanOperatorValue_spot = operatorValues(j,:)-meanRandSeqValue;
%calculate ka value from averaged curve
    y = meanOperatorValue_spot(associationFrames);
    x = timeArray(associationFrames);
    pka = polyfit(x, y, 1);
     ka_Osym(1,j) = pka(1);
end
kaOsymvalid = rmoutliers(ka_Osym,'median');
kaOsym = mean(kaOsymvalid);
kaOsym=kaOsym/lacIConcentration;





maxMutationNum = length(operatorSequence);
kaArray = nan(1,numel(uDnaSequences));
kdArray = nan(1,numel(uDnaSequences));
dnaArray = cell(1,numel(uDnaSequences));
kaCVArray = nan(1,numel(uDnaSequences));
kdCVArray = nan(1,numel(uDnaSequences));

for i=1:numel(uDnaSequences)
    if numel(uDnaSequences{i}) == maxMutationNum
        currMutationNum = sum(uDnaSequences{i}~=operatorSequence);
        if currMutationNum == mutationNum
            seqIndex = iu==i;
            operatorValues = cell2mat(spotIntensities(seqIndex));
%             meanOperatorValue = trimmean(operatorValues,30,1);
%             meanOperatorValue = meanOperatorValue-meanRandSeqValue;
            %calculate ka value from averaged curve
            ka_aDNA=nan(1,length(operatorValues(:,1)));
            kd_aDNA=nan(1,length(operatorValues(:,1)));
            for j = 1:length(operatorValues(:,1))
                meanOperatorValue_spot = operatorValues(j,:)-meanRandSeqValue;
            %calculate ka value from averaged curve
                y = meanOperatorValue_spot(associationFrames);
                x = timeArray(associationFrames);
                pka = polyfit(x, y, 1);
                y = meanOperatorValue_spot(disassociationFrames)/meanOperatorValue_spot(disassociationFrames(1));
                x = timeArray(disassociationFrames);
                pkd = polyfit(x, y, 1);
                ka_aDNA(1,j) = pka(1);
                kd_aDNA(1,j) = -pkd(1);
            end
            ka_valid = rmoutliers(ka_aDNA,'median')/lacIConcentration/kaOsym;
            kd_valid = rmoutliers(kd_aDNA,'median');
            kaStd = std(ka_valid);
            kdStd = std(kd_valid);
            ka_mean = mean(ka_valid);
            kd_mean = mean(kd_valid);
            ka_stderror = kaStd/sqrt(length(ka_valid));
            kd_stderror = kdStd/sqrt(length(kd_valid));
             if  ka_mean>=0 && kd_mean>=0 %&& ka_cv<1 && kd_cv<1 
                kaArray(i) = ka_mean;
                kdArray(i) = kd_mean;                
                dnaArray{i} = uDnaSequences{i};
                kaCVArray(i)= ka_stderror;
                kdCVArray(i)= kd_stderror;
            end
        end
    end
end
