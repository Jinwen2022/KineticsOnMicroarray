function [dataWildType,dataMutant]=sortByDNA(unsortedDataWildType,unsortedDataMutant)
[~,iWildType,iMutant] = intersect(unsortedDataWildType.dna,unsortedDataMutant.dna);
dataWildType.dna = unsortedDataWildType.dna(iWildType);
dataWildType.eq = unsortedDataWildType.eq(iWildType');
dataWildType.eqStderror = unsortedDataWildType.eqStderror(iWildType');
dataWildType.ka = unsortedDataWildType.ka(iWildType');
dataWildType.kaStderror = unsortedDataWildType.kaStderror(iWildType');
dataWildType.kd = unsortedDataWildType.kd(iWildType');
dataWildType.kdStderror = unsortedDataWildType.kdStderror(iWildType');
for i= 1:numel(unsortedDataWildType.mutTable(:,1,1))
    dataWildType.mutTable(i,:,:) = unsortedDataWildType.mutTable(iWildType(i),:,:);
end

dataMutant.dna = unsortedDataMutant.dna(iMutant);
dataMutant.eq = unsortedDataMutant.eq(iMutant');
dataMutant.eqStderror = unsortedDataMutant.eqStderror(iMutant');
dataMutant.ka = unsortedDataMutant.ka(iMutant');
dataMutant.kaStderror = unsortedDataMutant.kaStderror(iMutant');
dataMutant.kd = unsortedDataMutant.kd(iMutant');
dataMutant.kdStderror = unsortedDataMutant.kdStderror(iMutant');
for i = 1:numel(iMutant)
    dataMutant.mutTable(i,:,:) = unsortedDataMutant.mutTable(iMutant(i),:,:);
end

tf = strcmp(dataWildType.dna,dataMutant.dna);
if sum(tf) == numel(dataMutant.dna)
    fprintf('After intersection, both dataset were arranged by same DNA order\n')
else
    fprintf('Failed to arrange dataset by same DNA order')
end