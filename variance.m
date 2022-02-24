function varMeasure = variance(data,perc)

sorteddata=sort(data,'ascend');
numV=numel(sorteddata);
numToExcludeEachSide=round((perc/2)*numV);
startind=numToExcludeEachSide+1;
endind=numV-numToExcludeEachSide;

minval=sorteddata(startind);
maxval=sorteddata(endind);
meanCenter=mean(sorteddata(startind:endind));

varMeasure=(maxval-minval)/meanCenter;