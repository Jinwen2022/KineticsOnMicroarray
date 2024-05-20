function [tileCoord, tiles] = convertImageToTileCoord(imCoord,stitchROI)

stitchROI2 = reshape(stitchROI,[],size(stitchROI,3));
xmax = max(max(stitchROI2(:,2)));
ymax = max(max(stitchROI2(:,4)));
TL = zeros(ymax,xmax);
for i=1:size(stitchROI2,1)
    roi = stitchROI2(i,1:4);
    TL(roi(3):roi(4),roi(1):roi(2))=i;
end
tileCoord = nan(size(imCoord));
tiles = zeros(size(imCoord,1),1);
roundImCoord = round(imCoord);
for i=1:size(roundImCoord,1)
    s = roundImCoord(i,1:2);
    tiles(i) = TL(s(2),s(1));
    tileCoord(i,1:2)=s-stitchROI2(tiles(i),[1 3])+stitchROI2(tiles(i),[5 7]);
end
tileCoord(:,3:end)=imCoord(:,3:end);