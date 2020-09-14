function [orderedList, reachDistList, coreDistList, procesList] = faster_optics(data, Nmin, eps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for optics clustering of input 2d data
%improved version over the normal optics function

%INPUTS
% - data (2d data in the form of Nx2 where N is the number of data points)
% - Nmin (minimal number of points to form a cluster)
% - eps (maximum distance to check the Nmin parameter)

%OUTPUTS
% - orderedList (output indexing of the input data)
% - reachDistList (final reachability distance for every data point)
% - coreDistList (core distance for every data point)
% - procesList (if every data point was processed should be array of ones)

% to get the final reachabilty distance array with correct indexing use code simmilar to following 
%  for i = 1:size(data,1)
%     orderedReachList(i,1) = reachDistList(orderedList(i));
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reachDistList = zeros(size(data,1),1)-1; % -1 undefined
coreDistList = zeros(size(data,1),1)-1;
procesList = zeros(size(data,1),1);      % 0 unprocessed

orderedList = zeros(size(data,1),1);
orderedCount = 0;

seeds = zeros(size(data,1),1);

% neighIndices = zeros(size(data,1));
% neighCount = 0;

dist_matrix = pdist2(data,data);
unprocesCount = size(data,1);
seedCount_prev = 0;

while(unprocesCount)
    
    % get next unprocessed point
    i = getUnproces ( procesList ); 
    
    %get neighbors
    [neighIndices, neighCount] = getNeighbors ( i, eps, dist_matrix );
    
    % mark point as processed
    procesList(i) = 1;
    % increase order count, add point to ordered list
    orderedCount = orderedCount + 1;
    orderedList(orderedCount) = i;
    unprocesCount = unprocesCount - 1;
    
    % get core distance
    coreDistList(i) =  calcCoreDist ( neighIndices, i, neighCount, Nmin, coreDistList, dist_matrix);
    
    if coreDistList(i) ~= -1
        seeds(:) = -1;
        seedCount = 0;
        [seeds, seedCount, reachDistList, change] = update( neighIndices, neighCount, seeds, seedCount, i, procesList, reachDistList, coreDistList, dist_matrix);
        
        while ( seeds(1) ~= -1 )
            
            if seedCount_prev ~= seedCount || change

                seed_reach = zeros(seedCount,2);            

                seed_reach (:,1) = seeds(1:seedCount);
                seed_reach (:,2) = reachDistList(seeds(1:seedCount));

                seed_reach =  sortrows (seed_reach,2);
                seeds(1:seedCount) = seed_reach(:,1);
            end
            
            neighbor = seeds(1);
            
            procesList(neighbor) = 1;
            orderedCount = orderedCount + 1;
            orderedList(orderedCount) = neighbor;
            unprocesCount = unprocesCount - 1;
   
            seeds = seeds(2:end);
            seeds(end+1) = -1;
            seedCount = seedCount-1;
            
            [neighIndices_mark, neighCount_mark] = getNeighbors ( neighbor, eps, dist_matrix );
            
            coreDistList(neighbor) = calcCoreDist ( neighIndices_mark, neighbor, neighCount_mark, Nmin, coreDistList, dist_matrix);
            
            seedCount_prev = seedCount;
            
            if coreDistList(neighbor) ~= -1
                [seeds, seedCount, reachDistList, change] = update( neighIndices_mark, neighCount_mark, seeds, seedCount, neighbor, procesList, reachDistList, coreDistList, dist_matrix);
            end
        end
    end
    
            
            
            
            
end
end



function [indices, k] = getNeighbors ( i, eps, dist_matrix )
    % function that returns indices of all neighbors and their number (k)

    indices = find(dist_matrix(i,:)<=eps & dist_matrix(i,:)~=0);
    k = numel(indices);

end

function [coreDist] = calcCoreDist ( neighIndices, i, neighCount, Nmin, coreDistList, dist_matrix)
    % function that returns core distance for the selected point (i
    
    % if core distance was already calculated
    if coreDistList(i) ~= -1            % -1 represents undefined
        coreDist = coreDistList(i);
       
    else
    
    % if there is a minimal number of neighbors
        if neighCount >= Nmin-1
            distances = dist_matrix(i,neighIndices);
            distances = sort(distances);
            coreDist = distances(Nmin-1);
        else
        coreDist = -1;
        end
    end
end
   

function [seeds, seedNum, reachDistList, change] = update( neighIndices, neighCount, seeds, seedNum, i, procesList, reachDistList, coreDistList, dist_matrix)
%Update the seeds' reachability distance if a smaller value is found.

change = 0;
for k = 1:neighCount
    ind_neigh = neighIndices(k);
    
    
    if(procesList(ind_neigh)==0) %0 unprocessed
        new_reach = max ( [coreDistList(i); dist_matrix(i,ind_neigh)]);
        
        if(reachDistList(ind_neigh)==-1)
            reachDistList(ind_neigh) = new_reach;
            
            seedNum = seedNum +1;
            seeds(seedNum) = ind_neigh;
            
            
        else
            if new_reach < reachDistList(ind_neigh)
                reachDistList(ind_neigh) = new_reach;
                change = 1;
            end
        end
    end
end
end
      
function [indUnproc ] = getUnproces ( procesList )
% returns index of next unprocessed point

indUnproc = -1;
for k = 1: size(procesList,1)
    if procesList(k) == 0
        indUnproc = k;
        break
    end
end
end


    
    