function [ptots,koffmicros,kamax] = computePtotAndKoffmicro(kavals,kdvals,kamaxfrac,okinds)
if nargin<4
    okinds = true(size(kavals));
end

[sortedkavals,sortinds]=sort(kavals(okinds),'descend');
tempkdvals=kdvals(okinds);
sortedkdvals=tempkdvals(sortinds);
kavalscurr=sortedkavals(1:round(sum(okinds)*kamaxfrac));
kdvalscurr=sortedkdvals(1:numel(kavalscurr));

normX=std(kdvalscurr);
normY=std(kavalscurr);
relX=kdvalscurr/normX;
relY=kavalscurr/normY;
rng(1);
bestParams=fitLineErrorxy(relX,relY,1,1,1,1e-5,1e-5);
pfitCurr = [bestParams(1)*normY/normX bestParams(2)*normY];

kamax=pfitCurr(2);
ptots=kavals/kamax;
koffmicros=(kamax*kdvals)./(kamax-kavals);