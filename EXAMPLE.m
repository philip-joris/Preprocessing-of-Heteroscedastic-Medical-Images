% This script will illustrate how to set up and use the code in this
% repository. For an in-depth explanation of the algorithms and code,
% please see the original research article named:
%
%   "Preprocessing of Heteroscedastic Medical Images"
%
% Note: the provided .dcm image was obtained from the BrainWeb repository,
% which can be found at: http://brainweb.bic.mni.mcgill.ca/brainweb/
%
% Philip Joris
% 1/5/2018
% philip.joris@esat.kuleuven.be

%% Please first compile the cpp files using mex.

mex estimateISDs.cpp;
mex ISDClstr.cpp;

%% Load some data

% LOAD IMAGE (Change the path to example CT image)
img = double(squeeze(dicomread('brainweb_phantom.dcm')));

% Gaussian smoothing can be used to reduce noise or 'interpolate' discrete
% intensity images
%img = imgaussfilt3(img);

%% Estimate Intensity-specific distributions

nrOfBins        = 50;
[ISDs,edges]    = estimateISDs(img,nrOfBins);

%% EXAMPLE 1: Intensity transformation
% Note: CT images are far more heteroscedastic than MR images, so the
% effect of the transformation is more explicit when performed on CT
% images.

% CALCULATE INTENSITY-SPECIFIC STANDARD DEVIATIONS
centers = edges + (edges(2)-edges(1))/2;
centers = centers(1:end-1);
sigmas  = arrayfun(@(k) sqrt(sum((ISDs{k}-centers(k)).^2)/length(ISDs{k})), 1:nrOfBins);

% CALCULATE TRANSFORMATION
transform.x = edges;
transform.y = [0; cumsum(ones(size(sigmas(:)))./sigmas(:))];
transform.y = transform.y/max(transform.y);    

% TRANSFORM IMAGE (use interpolation)
imgT = reshape(interp1(transform.x,transform.y,img(:)),size(img));

figure; montage(permute(single(img),[1 2 4 3]),'DisplayRange',[]);
figure; montage(permute(single(imgT),[1 2 4 3]),'DisplayRange',[]);


%% EXAMPLE 2: Intensity classification
% Note: images are classified as a whole now. One can also first mask the
% image to a particular region (e.g. the brain region), and then cluster
% the remaining tissues.

% Calculate distance matrix
centers = edges + (edges(2)-edges(1))/2;
centers = centers(1:end-1);
sigmas  = arrayfun(@(k) sqrt(sum((ISDs{k}-centers(k)).^2)/length(ISDs{k})), 1:nrOfBins);
distanceMatrix = zeros(nrOfBins,nrOfBins);
for r=1:nrOfBins
    for c=1:nrOfBins
        distanceMatrix(r,c) = (centers(r)-centers(c)).^2/(sigmas(r)*sigmas(c));
    end
end

% Discretize image
imgD = discretize(img, edges);

% Cluster image
imgC = ISDClstr(int32(imgD), distanceMatrix, 5);

figure; montage(permute(single(imgC),[1 2 4 3]),'DisplayRange',[0 4]);
