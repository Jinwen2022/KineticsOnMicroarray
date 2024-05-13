function DataTable = writeDataLacI(mutTable,operatorSequences,legendStrs,ka,kd,eq,dna,kaStderror,kdStderror,eqStderror,nSpots)
DNAlabel = cell(length(ka),1);
labels_1mut = {'0_{sym}_1mut','0_1_1mut','0_2_1mut'};
labels_2mut = {'0_{sym}_2mut','0_1_2mut','0_2_2mut'};
    for i=1:numel(operatorSequences)
    ind =  mutTable(:,i,1);
    ind_1mut = mutTable(:,i,2);
    ind_2mut = mutTable(:,i,3);
    DNAlabel(ind)=legendStrs(i);
    DNAlabel(ind_1mut)=labels_1mut(i);
    DNAlabel(ind_2mut)=labels_2mut(i);
    end
DataTable = table(DNAlabel,ka',kd',eq',dna,kaStderror',kdStderror',eqStderror',nSpots','VariableNames',{'DNAtype','ka','kd','eq','dna','kaStderror','kdStderror','eqStderror','nSpots'});
end