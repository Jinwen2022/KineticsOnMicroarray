function [ka, p] = computeKa(data,time,frames,lacIConcentration,varargin)
ip = inputParser();
ip.addParameter('MakePlot',false);
ip.addParameter('Axes',[]);
ip.parse(varargin{:});
makePlot=ip.Results.MakePlot;
ax = ip.Results.Axes;

if nargin<5
    makePlot = false;
end

y = data(frames);
x = time(frames);
% if sum(y<0)>numel(y)/2
%     ka = nan;
%     return
% end

p = polyfit(x, y, 1);
ka = p(1)/lacIConcentration;

if makePlot
    if isempty(ax)
        figure, ax=gca;
    end
    f=polyval(p,x);
    plot(ax,x,y,'o',x,f,'-');
%     legend('data','linear fit')
end