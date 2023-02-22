function [ka, p] = computeKa_Cf(varargin)

ip = inputParser();
ip.addRequired('meanDNAValues',@(x) isnumeric(x));
ip.addRequired('spotIntensities',@(x) isnumeric(x));
ip.addRequired('time',@(x) isnumeric(x));
ip.addRequired('framesA',@(x) isnumeric(x));
ip.addRequired('lacIConcentration',@(x) isnumeric(x));
ip.addRequired('numC',@(x) isnumeric(x));
ip.parse(varargin{:});

meanDNAValues = ip.Results.meanDNAValues;
spotIntensities=ip.Results.spotIntensities;
time=ip.Results.time;
framesA=ip.Results.framesA;
lacIConcentration = ip.Results.lacIConcentration;
numC=ip.Results.numC;

y = spotIntensities(framesA);
x = time(framesA);
p = polyfit(x, y, 1);
ka = p(1)*numC/meanDNAValues/lacIConcentration/100;

end