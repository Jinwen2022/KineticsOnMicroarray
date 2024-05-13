function [mask,topLeftCorner, spotCentroids] = arrayGridMask(im,pxSize,topLeftCorner,debugFlag,contrast)
% Creates a binary mask of microarray spots for a microarray with either
% 96x164 or 266x170 spots.
% The spot radius is only 90% of the real radius in order to deal with
% the potential misalignment of the spot mask with images.
%
% Input:
%   im - preprocessed fluo image.
%   pxSize - pixel size in mum.
%   topLeftCorner - [x,y] coordinates of the top left spot centroid.
%                   Optional, if not provided or empty, the centroid is
%                   set using gui.  Default [].
%   debugFlag - optional, if true plot image with spot outlines. Default
%               false.

colSpacing = 0.127*1000/pxSize;
rowSpacing = 0.073323*1000/pxSize;
radius = round(0.9*10*pxSize);
if size(im,1)>rowSpacing*265 && size(im,2)>colSpacing*84
    NROWS = 266;
    NCOLS = 85;
elseif size(im,1)>rowSpacing*95 && size(im,2)>colSpacing*81
    NROWS = 96;
    NCOLS = 82;
else
    error('The image size is too small')
end

if nargin<3 || isempty(topLeftCorner)
    figure, imshow(im(1:round(10*rowSpacing),1:round(10*colSpacing)),[0,660]);
    [x,y] = ginput(1);
    topLeftCorner = [x y];
end
if nargin<4
    debugFlag = false;
end

spotCentroidsX = topLeftCorner(1):colSpacing:colSpacing*NCOLS;
spotCentroidsY = topLeftCorner(2):rowSpacing:rowSpacing*NROWS;
spotCentroids = zeros(NROWS,NCOLS*2,2);

mask = false(size(im));
[XX, YY] = meshgrid(spotCentroidsX, spotCentroidsY);
spotCentroids(:,1:2:end,1) = round(XX);
spotCentroids(:,1:2:end,2) = round(YY);
rcenters = round([XX(:) YY(:)]);
[XX, YY] = meshgrid(spotCentroidsX+colSpacing/2, spotCentroidsY+rowSpacing/2);
rcenters = [rcenters; round([XX(:) YY(:)])];
mask(sub2ind(size(mask),rcenters(:,2),rcenters(:,1))) = 1;
spotCentroids(:,2:2:end,1) = round(XX);
spotCentroids(:,2:2:end,2) = round(YY);
se = strel('disk',radius,8);
mask = imdilate(mask,se);
if debugFlag
<<<<<<< HEAD
    figure, imshow(im,[0,660]), hold on
=======
<<<<<<< HEAD
    figure, imshow(im,contrast), hold on
=======
    figure, imshow(im,[0,660]), hold on
>>>>>>> ec85aa6c343fe08d5df7c9a42a4c1d92269b8576
>>>>>>> 824d11ea3e0a815643ac501cfcfcf2f9f3572d60
    B = bwboundaries(mask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
end