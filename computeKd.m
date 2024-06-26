function [kd, fitObj] = computeKd(data,time,frames,varargin)
% Parse parameters
ip = inputParser;
ip.addOptional('Fittype','linear',@(x) ischar(x) && ismember(x,{'linear','exp'}));
ip.addParameter('StartPoint',[],@(x) isempty(x) || isnumeric(x))
ip.addParameter('Threshold',-Inf,@(x) isnumeric(x) && numel(x)==1)
ip.addParameter('MakePlot',false);
ip.addParameter('Axes',[]);
ip.parse(varargin{:});
threshold=ip.Results.Threshold;
fittype=ip.Results.Fittype;
makePlot=ip.Results.MakePlot;
ax = ip.Results.Axes;

if data(frames(1))<threshold
    kd = nan;
    return;
end

% Prepare data
y = data(frames)/data(frames(1));
x = time(frames)-time(frames(1));

% Do the fit
if strcmp(fittype,'linear')
    fitObj = polyfit(x, y, 1);
    kd = -fitObj(1);
elseif strcmp(fittype,'exp')
    try
        kd0 = ip.Results.StartPoint(1);
        if numel(ip.Results.StartPoint)==2
            a0 = ip.Results.StartPoint(2);
        else
            a0 = 1;
        end
        fitObj = fit(x(:),y(:), 'exp1','StartPoint',[a0 -kd0]);
        kd = -fitObj.b;
    catch
        kd = nan;
    end
else
    error(['Fit type ' fittype ' is not supported']); 
end

% Make plot if requested
if makePlot && isfinite(kd)
    if isempty(ax)
        figure, ax=gca;
    end
    if strcmp(fittype,'linear')
        f=polyval(fitObj,x);
    elseif strcmp(fittype,'exp')
        f=fitObj(x);
    end
    f=polyval(fitObj,x);
    plot(ax,x,y,'o',x,f,'-');
    hold on
    f=fitObj(x);
    plot(ax,x,y,'o',x,f,'-');
    if strcmp(fittype,'exp')
        hold(ax,'on')
        f0=exp(-kd0*x);
        plot(ax,x,f0,'b--');
        hold(ax,'off')
    end
end