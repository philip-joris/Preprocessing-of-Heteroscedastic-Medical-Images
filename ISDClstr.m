% -------------------------------------------------------------------------
% ISDClstr.m Help file for MEX-file
%
% Author:   Philip Joris
% Created:  1/5/2018
% Contact:  philip.joris@esat.kuleuven.be
%
% Clusters the intensites of the provided image into the specified number
% of groups. To cluster the intensities, the algorithm minimizes the within
% cluster distance (pairwise distance between all elements within a single
% cluster).
%
% The method follows the approach outlined in:
%
%   "Preprocessing of heteroscedastic medical images" 
%           by Joris. P et. al.
%
% Syntax of the method is:
%   labelImg = ISDClstr(image, distanceMatrix, k);
%
% where:
%
% INPUTS:
%
% > image:      3D numerical array of the type int32, where intensity
%               values correspond to bin labels. A intensity of value 1
%               corresponds to the first row/column in the distance matrix.
% > distanceMatrix:     2D numerical matrix of type double, where the entry
%                       at location (r,c) specifies the distance between an
%                       intensity in bin r and bin c. This distance matrix
%                       can be symmetrical or asymmetrical.
% > k:  Integer value specifying the number of clusters.
%
% OUTPUTS:
% 
% > labelImg:   3D numerical array of the type int32 and same size as the
%               input image. Each element in the array holds the label of
%               its corresponding cluster.
%
% -------------------------------------------------------------------------
