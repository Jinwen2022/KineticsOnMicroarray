function [ptots,koffmicros,kamax] = computePtotAndKoffmicro(kavals,kdvals,kamaxfrac,okinds)
if nargin<4
    okinds = true(size(kavals));
end

[sortedkavals,sortinds]=sort(kavals(okinds),'descend');
tempkdvals=kdvals(okinds);
sortedkdvals=tempkdvals(sortinds);
<<<<<<< HEAD
kavalscurr=sortedkavals(round(0.03*numel(sortedkavals)):round(sum(okinds)*kamaxfrac));%removed top 1% ka values that might be abnormal 
kdvalscurr=sortedkdvals(round(0.03*numel(sortedkavals)):round(sum(okinds)*kamaxfrac));
=======
kavalscurr=sortedkavals(1:round(sum(okinds)*kamaxfrac));%kon,max was estimated as the ka-intercept of a linear regression of the ka versus kd data, where the operators with the 5 % highest ka values were included in the regression
kdvalscurr=sortedkdvals(1:numel(kavalscurr));
>>>>>>> 8dfb5a1631a80d99bcb3dae40143bcc7a9849c16

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