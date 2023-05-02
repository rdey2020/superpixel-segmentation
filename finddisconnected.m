% Usage: list = finddisconnected(l)
%
% Argument:   l - A labeled image segmenting an image into regions, such as
%                 might be produced by a graph cut or superpixel algorithm.
%                 All pixels in each region are labeled by an integer.
%
% Returns: list - A cell array of lists of regions that are not
%                 connected. Typically there are 5 to 6 lists.
%

function list = finddisconnected(l)
 
    debug = 0;
    [Am, Al] = regionadjacency(l);
    
    N = max(l(:));  % number of labels
    
    % Array for keeping track of visited labels
    visited = zeros(N,1);

    list = {};
    listNo = 0;
    for n = 1:N

        if ~visited(n)
            listNo = listNo + 1;
            list{listNo} = n;
            visited(n) = 1;
            
            % Find all regions not directly connected to n and not visited
            notConnected = setdiff(find(~Am(n,:)), find(visited));
            
            % For each unconnected region check that it is not already
            % connected to a region in the list. If not, add to list
            for m = notConnected
                if isempty(intersect(Al{m}, list{listNo}))
                    list{listNo} = [list{listNo} m];
                    visited(m) = 1;
                end
            end
         end % if not visited(n)
        
    end
    
    % Display each list of unconncted regions as an image
    if debug   
        for n = 1:length(list)
            
            mask = zeros(size(l));
            for m = 1:length(list{n})
                mask = mask | l == list{n}(m);
            end
            
            fprintf('list %d of %d length %d \n', n, length(list), length(list{n}))
            show(mask);
            keypause
        end
    end
    