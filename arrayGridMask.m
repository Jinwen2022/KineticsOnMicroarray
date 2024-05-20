function [mask,topLeftCorner, spotCentroids, bkgMask] = arrayGridMask(im,pxSize,topLeftCorner,debugFlag,contrast)
% Creates a labelled mask of microarray spots for a microarray with either
% 96x164 or 266x170 spots, and optionally a labelled mask for local
% background subtraction.
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
%   contrast - a tuple with upper and lower pixel intensity value cut-offs
%              for image plotting. By default, [0 660].
%
% Output:
%   mask - an image of the same size as im with segmented and labelled
%          microarray spots;
%   topLeftCorner - (x,y)-coordinates of the top left microarray spot
%                   centroid;
%   spotCentroids - coordinates of all microarray spot centroids
%   bkgMask - optional, a labelled mask of background pixels around
%             microarray spots


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
if nargin<5
    contrast = [0 660];
end

spotCentroidsX = topLeftCorner(1):colSpacing:colSpacing*NCOLS;
spotCentroidsY = topLeftCorner(2):rowSpacing:rowSpacing*NROWS;
spotCentroids = zeros(NROWS,NCOLS*2,2);

ceedMask = zeros(size(im));
[XX, YY] = meshgrid(spotCentroidsX, spotCentroidsY);
spotCentroids(:,1:2:end,1) = XX;
spotCentroids(:,1:2:end,2) = YY;
rcenters = round([XX(:) YY(:)]);
[XX, YY] = meshgrid(spotCentroidsX+colSpacing/2, spotCentroidsY+rowSpacing/2);
spotCentroids(:,2:2:end,1) = XX;
spotCentroids(:,2:2:end,2) = YY;
rcenters = [rcenters; round([XX(:) YY(:)])];
[~,ind] = sort(rcenters(:,1));
rcenters = rcenters(ind,:);
ceedMask(sub2ind(size(ceedMask),rcenters(:,2),rcenters(:,1))) = 1:numel(spotCentroids)/2;

se = strel('disk',radius,8);
mask = imdilate(ceedMask,se);
if debugFlag
    figure, imshow(im,contrast), hold on
    B = bwboundaries(mask,'noholes');
    visboundaries(B,'EnhanceVisibility',0,'LineWidth',1);
end

if nargout==4
    diskSe = strel('disk',round(1.35*10*pxSize),8);
    outMask = imdilate(ceedMask,diskSe);
    rectSe = strel('rectangle',odd([round(rowSpacing*0.95) size(diskSe.Neighborhood,1)+3]));
    bkgMask = imdilate(ceedMask,rectSe);
    bkgMask = bkgMask-outMask;
    if debugFlag
        figure, imshow(labeloverlay(imadjust(im),bkgMask)),axis image
    end
end