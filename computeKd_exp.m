function [kd_exp] = computeKd_exp(data,time,frames,varargin)
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

% Prepare data for fit
% Do the linear fit
y = data(frames)/data(frames(1));
x = time(frames)-time(frames(1));
%get intial slope for linear fit
y_linear = y(1:8);
x_linear = x(1:8);
fitObj_linear = polyfit(x_linear, y_linear, 1);
kd_linear = -fitObj_linear(1);
%Do the exp fit
    try
        % use the initial slope approximation for exponential fit starting
        % point
        kd0 = kd_linear;
        a0 = kd_linear*2;
        fitObj_exp = fit(x(:),y(:), 'exp1','StartPoint',[a0 -kd0]);
        kd_exp = -fitObj_exp.b;
    catch
        kd_exp = nan;
    end
    


% Make plot if requested
if makePlot && isfinite(kd_linear)
    if isempty(ax)
        figure, ax=gca;
    end
    
    f_linear=polyval(fitObj_linear,x_linear);

    f_exp=fitObj_exp(x);
    
    plot(ax,x,y,'o',x_linear,f_linear,'b-');
    hold(ax,'on')
    plot(ax,x,f_exp,'r-');
    
%     hold(ax,'on')
%     f0=exp(-kd0*x);
%     plot(ax,x,f0,'r--');
    hold(ax,'off')
    legend('data','linear fit','exp fit','exp starting point with kd0 from linear fit')
    txt = {'Linear fit kd:' kd_linear,'exponential fit kd:' kd_exp};
    text(2000,0.7,txt)
end
end