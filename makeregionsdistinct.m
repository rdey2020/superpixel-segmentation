% Usage: [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
%
% Arguments: seg - A region segmented image, such as might be produced by a
%                  superpixel or graph cut algorithm.  All pixels in each
%                  region are labeled by an integer.
%   connectivity - Optional parameter indicating whether 4 or 8 connectedness
%                  should be used.  Defaults to 4.ͨ
%
% Returns:   seg - A labeled image where all segments are distinct.
%       maxlabel - Maximum segment label number.  
%

function [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
    
    if ~exist('connectivity', 'var'), connectivity = 4; end
    
    % Ensure every segment is distinct but do not touch segments 
    % with a label of 0
    labels = unique(seg(:))'; 
    maxlabel = max(labels);
    labels = setdiff(labels,0);  
    
    for l = labels
        [bl,num] = bwlabel(seg==l, connectivity);    
        if num > 1  % We have more than one region with the same label
            for n = 2:num
                maxlabel = maxlabel+1;  % Generate a new label
                seg(bl==n) = maxlabel;  % and assign to this segment
            end
        end
    end