function [im,angularDisplacement,cropRoi] = preprocessImage(im,angularDisplacement,cropRoi)
% Performs microarray image rotation and cropping.
%
% Input:
%   im - microarray image.
%   angularDisplacement - scalar, defines counter-clockwise image rotation
%                         in degrees. If empty, a gui is used to define
%                         angularDisplacement. Default [].
%   cropRoi - either a 4-element vector [XMIN YMIN WIDTH HEIGHT] for image
%             cropping region, or empty. If empty, a gui is used to define
%             a cropping region. Default [].

if nargin<2
    angularDisplacement = [];
end
if nargin<3
    cropRoi = [];
end

if isempty(angularDisplacement)
    angularDisplacement = getAngularDisplacement(imadjust(im));
end
im = imrotate(im,angularDisplacement,'bilinear');

if isempty(cropRoi)
    cropRoi = getROI(imadjust(im),'Crop array');
end
im = imcrop(im,cropRoi);