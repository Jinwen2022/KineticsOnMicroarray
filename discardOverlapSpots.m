function correctedDnaSequences = discardOverlapSpots(dnaSequences,dirCy3,parameters)

[imRaw, stitchROI] = stitchArrayImages(dirCy3{1},parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,1,-1);
im = preprocessImage(imRaw,parameters.angularDisplacement,parameters.roi);
[~,~,Trot] = getAngularDisplacement(imRaw,parameters.angularDisplacement);
Tcrop = getTranslationMatrix(-parameters.roi+1);
T = Trot * Tcrop;
[~,~,centroids] = arrayGridMask(im,parameters.pxSize,parameters.topLeftCorner,0);
centroidsRaw = convertCoord(reshape(centroids,[],2),T,1);

discardSpots = false(size(dnaSequences));
bordersX= unique(stitchROI(:,:,1));
bordersY= unique(stitchROI(:,:,3));
spotRadius = 15;
for i=1:size(centroidsRaw,1)
    currSpot = centroidsRaw(i,:);
    if any(currSpot(1)>bordersX-spotRadius & currSpot(1)<bordersX+spotRadius) || ...
       any(currSpot(2)>bordersY-spotRadius & currSpot(2)<bordersY+spotRadius)
        discardSpots(i)= true;
    end
end
correctedDnaSequences = dnaSequences;
correctedDnaSequences(discardSpots)={''};