% Usage: [nL, minLabel, maxLabel] = renumberregions(L)
%
% Argument:   L - A labeled image segmenting an image into regions, such as
%                 might be produced by a graph cut or superpixel algorithm.
%                 All pixels in each region are labeled by an integer.
%
% Returns:   nL - A relabeled version of L so that label numbers form a
%                 sequence 1:maxRegions  or 0:maxRegions-1 depending on
%                 whether L has a region labeled with 0s or not.
%      minLabel - Minimum label in the renumbered image.  This will be 0 or 1.
%      maxLabel - Maximum label in the renumbered image.
%

function [nL, minLabel, maxLabel] = renumberregions(L)

    nL = L;
    labels = unique(L(:))';  % Sorted list of unique labels 
    N = length(labels);
    
    % If there is a label of 0 we ensure that we do not renumber that region
    % by removing it from the list of labels to be renumbered.
    if labels(1) == 0
        labels = labels(2:end);
        minLabel = 0;
        maxLabel = N-1;
    else
        minLabel = 1;
        maxLabel = N;
    end
    
    % Now do the relabelling
    count = 1;
    for n = labels
        nL(L==n) = count;
        count = count+1;
    end
    