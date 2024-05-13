function [myLines] = plotDissociationCurves(varargin)
% Similar to plotAverageSpotIntensity but plot only dissociation curves
% and require the starting frame of dissociation to normalize fluorescence
% curves by the fluorescence at the beginning of dissociation.
%% Parse parameters
if length(varargin)>1 && isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = varargin{1};
    varargin = varargin(2:end);
else
    figure
    ax = gca;
end
ip = inputParser();
ip.addRequired('dnaSequences',@(x) iscell(x));
ip.addRequired('spotIntensities',@(x) iscell(x));
ip.addRequired('operatorSequences',@(x) iscell(x) || ischar(x) || isstring(x));
ip.addRequired('randSequence',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addRequired('startFrame',@(x) isnumeric(x));
ip.addOptional('timeArray',[],@(x) isnumeric(x));
ip.addParameter('Color','', @(x) iscell(x) || ischar(x) || isnumeric(x));
ip.addParameter('LineStyle','-',@(x) ischar(x) | isstring(x));
ip.addParameter('LineWidth',1,@(x) isnumeric(x));
ip.addParameter('alphaScale',1,@(x) isnumeric(x) & (x>=0) & (x<=1));
ip.addParameter('Percent',30,@(x) isnumeric(x) && numel(x)==1 && x>=0 && x<=100);
ip.addParameter('Handle',false,@(x) islogical(x) | (x==0) | (x==1));

ip.parse(varargin{:});

dnaSequences = ip.Results.dnaSequences;
spotIntensities = ip.Results.spotIntensities;
operatorSequences = ip.Results.operatorSequences;
randSequence = ip.Results.randSequence;
startFrame = ip.Results.startFrame;
timeArray = ip.Results.timeArray;
colors = ip.Results.Color;
linestyle = ip.Results.LineStyle;
linewidth = ip.Results.LineWidth;
alphaScale = ip.Results.alphaScale;
percent = ip.Results.Percent;
linesHandle = ip.Results.Handle;

if ~iscell(operatorSequences)
    operatorSequences = {operatorSequences};
end
if isempty(timeArray)
    timeArray = 1:numel(spotIntensities{1});
end
if isempty(colors)
    colors = lines(numel(operatorSequences));
    colors = mat2cell(colors,ones(size(colors,1),1));
end

if linesHandle
    myLines=[];
end
%% Main code
if isempty(randSequence)
    meanRandSeqValues = 0;
else
    randSeqIndex = strcmp(dnaSequences,randSequence);
    randSeqArrayValues = cell2mat(spotIntensities(randSeqIndex));
    meanRandSeqValues = trimmean(randSeqArrayValues,percent,1);
end
hold on
dissocTime = timeArray(startFrame:end)-timeArray(startFrame);
for i=1:numel(operatorSequences)
    operatorIndex = strcmp(dnaSequences,operatorSequences{i});
    if any(any(operatorIndex))
        operatorArrayValues = cell2mat(spotIntensities(operatorIndex));
        meanOperatorValues = trimmean(operatorArrayValues,percent,1);
        meanOperatorValues = meanOperatorValues-meanRandSeqValues;
        dissocValues = meanOperatorValues(startFrame:end)/meanOperatorValues(startFrame);
        if linesHandle
             p = plot(ax,dissocTime,dissocValues,linestyle,'Color',colors{i},'LineWidth',linewidth);
             myLines(end+1) = p;
        else
            p=plot(ax,dissocTime,dissocValues,linestyle,'Color',colors{i},'LineWidth',linewidth);
        end
        if alphaScale<1
            p.Color(4) = alphaScale;
        end
    end
end