function [myLines,outMeans, outSems] = plotAverageSpotIntensity(varargin)
% Computes mean spot intensity per operator in each frame. Optionally,
% subtracts spot intensity of a given random sequence.
% The function plots resulting average spot intensity curves per operator
% and returns the mean values and the standard error of mean.
%
% Input:
%   dnaSequences - 2D cell array with DNA sequences in microarray spots.
%   spotIntensities - 2D cell array, each element is a vector with average
%                     spot intensities of the corresponding microarray spot
%                     over acquisition frames.
%   operatorSequences - cell array with operator sequences
%   randSequence - char array with a random DNA sequences used for signal
%                  subtraction. If '', skip subtraction.
%   timeArray - optional, vector with frame acqusition times which are used
%               for the x-axis in the plot. If missing, the frame numbers
%               are used instead.
%
% plotAverageSpotIntensity(AX,...) plots into the axes with handle AX.
%
% plotAverageSpotIntensity(...,'Errorbar',flag) uses standar error of mean
% as errorbars along y-axis if flag is true. Default false.
%
% plotAverageSpotIntensity(...,'Color',colors) sets operator-specific line
% colors defined by cell array 'colors'. By default, 'lines' colormap is
% used.
%
% plotAverageSpotIntensity(...,'LineStyle',linestyle) sets all line styles
% as defined by char array linestyle. Default '-'.
%
% plotAverageSpotIntensity(...,'LineWidth',linewidth) sets all line widths
% in points as specified by a positive value linewidth. Default 1.
%
% plotAverageSpotIntensity(...,'alphaScale',alphascale) controls line 
% opacity, alphascale must be between 0 and 1 (full opacity). Default 1.
%
% plotAverageSpotIntensity(...,'Percent',perc) is used to compute trimmed
% average spot intensity per operator and frame, see trimmean. Default 30.
%
%
% Output:
%   outMeans - cell array, were element i contains average spot
%              intensities of operator operatorSequences{i}
%   outSems  - cell array, were element i contains standard error of mean
%              of spot intensities of operator operatorSequences{i}.
%              Computed only if the errorbar flag is true.
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
ip.addRequired('spotIntensities',@(x) isnumeric(x));
ip.addRequired('operatorSequences',@(x) iscell(x) || ischar(x) || isstring(x));
ip.addRequired('randSequence',@(x) isempty(x) || ischar(x) || isstring(x));
ip.addOptional('timeArray',[],@(x) isnumeric(x));
ip.addParameter('Errorbar',false,@(x) islogical(x) | (x==0) | (x==1));
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
timeArray = ip.Results.timeArray;
colors = ip.Results.Color;
errorbarFlag = ip.Results.Errorbar;
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
    randSeqArrayValues = spotIntensities(:,randSeqIndex);
    meanRandSeqValues = trimmean(randSeqArrayValues,percent,1);
end
if errorbarFlag || nargout==3
    sems = cell(size(operatorSequences));
end
hold on
meanOperatorValues = cell(size(operatorSequences));
for i=1:numel(operatorSequences)
    operatorIndex = strcmp(dnaSequences,operatorSequences{i});
    if any(any(operatorIndex))
        operatorArrayValues = spotIntensities(:,operatorIndex);
        meanOperatorValues{i} = trimmean(operatorArrayValues,percent,2);
        meanOperatorValues{i} = meanOperatorValues{i}-meanRandSeqValues;
        if errorbarFlag || nargout==3
            sems{i} = std(operatorArrayValues,0,1)/sqrt(size(operatorArrayValues,2));
        end
        if errorbarFlag
            p=errorbar(ax,timeArray,meanOperatorValues{i},sems{i},sems{i},linestyle,'Color',colors{i},'LineWidth',linewidth);
        else
            if linesHandle
                p=plot(ax,timeArray,meanOperatorValues{i},linestyle,'Color',colors{i},'LineWidth',linewidth);
                myLines(end+1) = p;
            else
                p=plot(ax,timeArray,meanOperatorValues{i},linestyle,'Color',colors{i},'LineWidth',linewidth);
            end
        end
        if alphaScale<1
            p.Color(4) = alphaScale;
        end
    end
end
% hold off
%% Handle output if required
if nargout
    outMeans = meanOperatorValues;
    if nargout==3
        outSems = sems;
    end
end