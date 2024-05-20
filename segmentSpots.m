function [mask, ceedsRaw] = segmentSpots(im,thr,displayFlag,clim)

if nargin<2
    thr=1.5;
end
if nargin<3
    displayFlag=0;
end
if nargin<4
    clim=[];
end

% set parameters
gaussSigma = 4;
spotRadius = 15;
minSpotArea = 1000;
maxSpotArea = 20000;
minExtent = 0.5;
maxAxisRatio  = 1.5;

% compute magnitude
im0 = im;
im = imgaussfilt(im,gaussSigma);
[imgx, imgy] = BuildTrajectories.derivative5(im, 'x', 'y');
mag = sqrt(imgx.^2 + imgy.^2)+eps; % (+eps to avoid division by 0)

% threshold the magnitude image and do some 
bw = mag>thr;
bw = bw & ~bwareaopen(bw,maxSpotArea);
bw = imfill(bw,"holes");
bw = bwareaopen(bw,minSpotArea);
% split connected spots
bw = imerode(bw,strel('disk',round(0.8*spotRadius),4));

% find spot centroids and create spot mask
stats = regionprops(bw,im,'Extent','Centroid','MajorAxisLength','MinorAxisLength');
keepMask = ([stats.Extent]>minExtent) & ...
    ([stats.MajorAxisLength]./[stats.MinorAxisLength] < maxAxisRatio);
ceedsRaw = reshape([stats(keepMask).Centroid],2,[])';
ceeds = round(ceedsRaw);
mask = false(size(im));
mask(sub2ind(size(mask),round(ceeds(:,2)),round(ceeds(:,1))))=true;
mask = imdilate(mask,strel('disk',spotRadius,8));

if displayFlag
    figure
    if isempty(clim)
        imagesc(imadjust(im0));
        %imagesc(labeloverlay(imadjust(im0),mask,'Transparency',0.7));
    else
        imagesc(im0,clim)
        %imagesc(labeloverlay(im0,mask,'Transparency',0.7),clim);
    end
    colormap gray
    hold on
    visboundaries(bwboundaries(mask,'noholes'),'EnhanceVisibility',0,'Color','r')
    axis image
end