function plotKaVsKaKdvsKd(dna1,dna2,ka1,ka2,kd1,kd2,mutTable1,mutTable2,labels)
if nargin<9
    labels = {'','','',''};
end

ind = ~cellfun('isempty',labels);
labels(ind) = strcat({' '},labels(ind));
if numel(labels) == 2 
    labels = repmat(labels,size(labels));
end

[~,ind1,ind2] = intersect(dna1,dna2);

ax=nexttile;
kaOsym1 = ka1(mutTable1(:,1,1));
kaOsym2 = ka2(mutTable2(:,1,1));
scatter(ax,ka1(ind1)/kaOsym1,ka2(ind2)/kaOsym2);
xlim([0,1.2])
ylim([0,1.2])
pbaspect([1 1 1])
xlabel(['ka/ka(Osym)' labels{1}])
ylabel(['ka/ka(Osym)' labels{2}])

ax=nexttile;
scatter(ax,kd1(ind1),kd2(ind2));
xlim([0,1.5e-3])
ylim([0,1.5e-3])
pbaspect([1 1 1])
xlabel(['kd (s^{-1})' labels{3}])
ylabel(['kd (s^{-1})' labels{4}])
