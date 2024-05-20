function spotIntensities = computeSpotIntensities(im,mask,bkgMask)
% Compute average spot fluorescence insensities in microarray image *im* 
% using a binary image *mask* with spot locations. If *bkgMask* is provided
% do local background subtraction.
% Return either 96x164 or 266x170 matrix with spot intensities.
props = regionprops(mask,im,'MeanIntensity');
if numel(props) == 96*164
    nRows = 96;
    nCols = 164;
elseif numel(props) == 266*170
    nRows = 266;
    nCols = 170;
else
    error('There must be 96x164 or 266x170 microarray spots')
end
spotIntensities = [props.MeanIntensity];
if nargin==3 % Apply local background subtraction
    props = regionprops(bkgMask,im,'MeanIntensity');
    localBkg = [props.MeanIntensity];
    spotIntensities = spotIntensities - localBkg;
end
spotIntensities = reshape(spotIntensities,nRows,nCols);
