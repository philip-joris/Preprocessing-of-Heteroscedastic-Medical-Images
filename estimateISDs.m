% -------------------------------------------------------------------------
% estimateISDs.m Help file for MEX-file
%
% Author:   Philip Joris
% Created:  1/5/2018
% Contact:  philip.joris@esat.kuleuven.be
%
% Estimates the intensity-specific distributions for the provided image.
% Multiple images can be concatenated along the 3rd dimension to learn the
% ISDs from multiple images simultaneously (for CT images for example). The
% number of bins needs to be specified as well.
%
% The method follows the approach outlined in:
%
%   "Preprocessing of heteroscedastic medical images" 
%           by Joris. P et. al.
%
% Syntax of the method is:
%   [ISDs,edges] = estimateISDs(image, nrOfBins);
%
% where:
%
% INPUTS:
%
% > image:      3D numerical array of the type double.
% > nrOfBins:   Specifies the number of bins used to discretise the image.
%
% OUTPUTS:
% 
% > ISDs:   Cell array of length 'nrOfBins', where each cell element
%           contains the ISD for a specific bin (intensity).
% > edges:  Column vector containing the edges of the bins that were used
%           to discretize the image.
%
% -------------------------------------------------------------------------
