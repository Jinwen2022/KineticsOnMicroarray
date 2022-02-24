function spotIntensities = computeSpotIntensities(im,mask)
% Compute average spot fluorescence insensities in microarray image *im* 
% using a binary image *mask* with spot locations.
% Return either 96*164 or 266*170 matrix with spot intensities.
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
spotIntensities = reshape([props.MeanIntensity],nRows,nCols);