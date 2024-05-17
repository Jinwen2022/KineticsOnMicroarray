function [ptots,koffmicros,kamax] = computePtotAndKoffmicro(kavals,kdvals,kamaxfrac,okinds)
if nargin<4
    okinds = true(size(kavals));
end

[sortedkavals,sortinds]=sort(kavals(okinds),'descend');
tempkdvals=kdvals(okinds);
sortedkdvals=tempkdvals(sortinds);
kavalscurr=sortedkavals(round(0.03*numel(sortedkavals)):round(sum(okinds)*kamaxfrac));%removed top 1% ka values that might be abnormal 
kdvalscurr=sortedkdvals(round(0.03*numel(sortedkavals)):round(sum(okinds)*kamaxfrac));

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