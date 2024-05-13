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
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.01]));
=======
<<<<<<< HEAD
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.05]));
=======
<<<<<<< HEAD
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.05]));
=======
<<<<<<< HEAD
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.05]));
=======
    angularDisplacement = getAngularDisplacement(imadjust(im,[0,0.01]));
>>>>>>> dba439a7c8255a67193fff196b5b04e036733833
>>>>>>> 1922fda1d52536b4bf47535f496ab3c8369ffdf1
>>>>>>> ec85aa6c343fe08d5df7c9a42a4c1d92269b8576
>>>>>>> 824d11ea3e0a815643ac501cfcfcf2f9f3572d60
end
im = imrotate(im,angularDisplacement,'bilinear');

if isempty(cropRoi)
<<<<<<< HEAD
    cropRoi = getROI(imadjust(im,[0,2000/65535]),'Crop array');
=======
<<<<<<< HEAD
    cropRoi = getROI(imadjust(im,[0,0.05]),'Crop array');
=======
<<<<<<< HEAD
    cropRoi = getROI(imadjust(im,[0,0.05]),'Crop array');
=======
<<<<<<< HEAD
    cropRoi = getROI(imadjust(im,[0,0.05]),'Crop array');
=======
    cropRoi = getROI(imadjust(im,[0,0.01]),'Crop array');
>>>>>>> dba439a7c8255a67193fff196b5b04e036733833
>>>>>>> 1922fda1d52536b4bf47535f496ab3c8369ffdf1
>>>>>>> ec85aa6c343fe08d5df7c9a42a4c1d92269b8576
>>>>>>> 824d11ea3e0a815643ac501cfcfcf2f9f3572d60
end
im = imcrop(im,cropRoi);