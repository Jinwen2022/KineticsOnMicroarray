function imCoord = convertTileToImageCoord(tileCoord,stitchROI, tiles)

stitchROI2 = reshape(stitchROI,[],size(stitchROI,3));
xmax = max(max(stitchROI2(:,2)));
ymax = max(max(stitchROI2(:,4)));
TL = zeros(ymax,xmax);
for i=1:size(stitchROI2,1)
    roi = stitchROI2(i,1:4);
    TL(roi(3):roi(4),roi(1):roi(2))=i;
end
imCoord = nan(size(tileCoord));
roundTileCoord = round(tileCoord);
for i=1:size(roundTileCoord,1)
    s = roundTileCoord(i,:);
    imCoord(i,1:2)=s+stitchROI2(tiles(i),[1 3])-stitchROI2(tiles(i),[5 7]);
end