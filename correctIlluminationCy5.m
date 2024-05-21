function W = correctIlluminationCy5(imDir,parameters, plotFlag)

[imRaw, stitchROI] = stitchArrayImages(imDir,parameters.overlap,parameters.angularDisplacementTile,parameters.rangePositions,1,-1);
im = preprocessImage(imRaw,parameters.angularDisplacement,parameters.roi);
[~,~,Trot] = getAngularDisplacement(imRaw,parameters.angularDisplacement);
Tcrop = getTranslationMatrix(-parameters.roi+1);
T = Trot * Tcrop;

% compute background mask
[~, spots] = segmentSpots(im,1.5,0);
ceedMask = zeros(size(im));
rcenters = round(spots);
ceedMask(sub2ind(size(ceedMask),rcenters(:,2),rcenters(:,1))) = 1:size(spots,1);
diskSe = strel('disk',round(1.35*10*parameters.pxSize),8);
outMask = imdilate(ceedMask,diskSe);
rowSpacing = 0.073323*1000/parameters.pxSize;
rectSe = strel('rectangle',odd([round(rowSpacing*0.95) size(diskSe.Neighborhood,1)+3]));
bkgMask = imdilate(ceedMask,rectSe);
bkgMask = bkgMask-outMask;

% Cy5 weights
stats = regionprops(bkgMask,im,'MeanIntensity');
spotsCoordsRaw = convertCoord(spots,T,1);
fitdata = [spotsCoordsRaw [stats.MeanIntensity]'];
fitdata = convertImageToTileCoord(fitdata,stitchROI);
fitObj = fit(fitdata(:,1:2),fitdata(:,3),'poly21');
imAdapterObj = genericReadAsFrames('metadata.txt',imDir);
posList = imAdapterObj.getPositionList();
chanNames = imAdapterObj.getChannels();
imInfo = imfinfo(fullfile(imDir,posList{1},chanNames{1},'img_000000000.tiff'));

% subtract camera offset and compute weigths
[X,Y]=meshgrid(1:imInfo.Width,1:imInfo.Height);
cameraOffset = 100;
Z = fitObj(X,Y)-cameraOffset;
W = max(max(Z))./Z;

if plotFlag
    figure, imshow(labeloverlay(imadjust(im),bkgMask)),axis image
    drawnow

    figure, 
    scatter3(fitdata(:,1),fitdata(:,2),fitdata(:,3),'b.')
    xlabel('x'),ylabel('y'),zlabel('w')
    hold on
    [X,Y] = meshgrid(1:10:imInfo.Width,1:10:imInfo.Height);
    surf(X,Y,fitObj(X,Y))
    zlim(quantile(fitdata(:,3),[0 0.995]))
    drawnow
    
    figure,
    imagesc(W)
    axis image
    colorbar
end