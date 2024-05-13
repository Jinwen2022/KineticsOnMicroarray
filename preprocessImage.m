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
<<<<<<< HEAD
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.05]));
=======
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.01]));
>>>>>>> dba439a7c8255a67193fff196b5b04e036733833
end
im = imrotate(im,angularDisplacement,'bilinear');

if isempty(cropRoi)
<<<<<<< HEAD
    cropRoi = getROI(imadjust(im,[0,0.05]),'Crop array');
=======
    cropRoi = getROI(imadjust(im,[0,0.01]),'Crop array');
>>>>>>> dba439a7c8255a67193fff196b5b04e036733833
end
im = imcrop(im,cropRoi);