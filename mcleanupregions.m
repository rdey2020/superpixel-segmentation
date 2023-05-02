% Usage: [seg, Am] = mcleanupregions(seg, seRadius)
%
% Arguments: seg - A region segmented image, such as might be produced by a
%                  graph cut algorithm.  All pixels in each region are labeled
%                  by an integer. 
%       seRadius - Structuring element radius.  This can be set to 0 in which
%                  case  the function will simply ensure all labeled regions
%                  are distinct and relabel them if necessary. 

function [seg, Am, mask] = mcleanupregions(seg, seRadius)
option = 2;
    [seg, maxlabel] = makeregionsdistinct(seg);
    
    % Perform a morphological opening on each segment, subtract the opening
    % from the orignal segment to obtain regions to be reassigned to
    % neighbouring segments.
    if seRadius
        se = circularstruct(seRadius);   
        mask = zeros(size(seg));

        if option == 1        
            for l = 1:maxlabel
                b = seg == l;
                mask = mask | (b - imopen(b,se));
            end
            
        else   
            list = finddisconnected(seg);
            
            for n = 1:length(list)
                b = zeros(size(seg));
                for m = 1:length(list{n})
                    b = b | seg == list{n}(m);
                end

                mask = mask | (b - imopen(b,se)); 
            end
        end
        
        % Compute distance map on inverse of mask
        [~, idx] = bwdist(~mask);
        
        seg(mask) = seg(idx(mask));
    end
    
    seg = makeregionsdistinct(seg);
    [seg, minLabel, maxLabel] = renumberregions(seg);
    Am = regionadjacency(seg);    
    