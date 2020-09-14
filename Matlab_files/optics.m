function [orderedList, reachDistList, coreDistList, procesList] = optics(data, Nmin, eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for optics clustering of input 2d data

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

unprocesCount = size(data,1);

while(unprocesCount)
    
    % get next unprocessed point
    i = getUnproces ( procesList ); 
    
    %get neighbors
    [neighIndices, neighCount] = getNeighbors (data, i, eps );
    
    % mark point as processed
    procesList(i) = 1;
    % increase order count, add point to ordered list
    orderedCount = orderedCount + 1;
    orderedList(orderedCount) = i;
    unprocesCount = unprocesCount - 1;
    
    % get core distance
    coreDistList(i) =  calcCoreDist (data, neighIndices, i, neighCount, Nmin, coreDistList);
    
    if coreDistList(i) ~= -1
        seeds(:) = -1;
        seedCount = 0;
        [seeds, seedCount, reachDistList] = update(data, neighIndices, neighCount, seeds, seedCount, i, procesList, reachDistList, coreDistList);
        
        while ( seeds(1) ~= -1 )
            
%           seeds(seedCount+1) = -1;
            seed_reach = 0;            
            
            seed_reach = zeros(seedCount,2);            

            seed_reach (:,1) = seeds(1:seedCount);
            seed_reach (:,2) = reachDistList(seeds(1:seedCount));

            seed_reach =  sortrows (seed_reach,2);
            seeds(1:seedCount) = seed_reach(:,1);
            
            neighbor = seeds(1);
            
            procesList(neighbor) = 1;
            orderedCount = orderedCount + 1;
            orderedList(orderedCount) = neighbor;
            unprocesCount = unprocesCount - 1;
   
            seeds = seeds(2:end);
            seeds(end+1) = -1;
            seedCount = seedCount-1;
            
            [neighIndices_mark, neighCount_mark] = getNeighbors (data, neighbor, eps );
            
            coreDistList(neighbor) = calcCoreDist (data, neighIndices_mark, neighbor, neighCount_mark, Nmin, coreDistList);
            
            if coreDistList(neighbor) ~= -1
                [seeds, seedCount, reachDistList] = update(data, neighIndices_mark, neighCount_mark, seeds, seedCount, neighbor, procesList, reachDistList, coreDistList);
            end
        end
    end
    
            
            
            
            
end
end


function dist = euclidianDist( x1, x2, y1, y2 )
    % function that returns euclidian distance between two points
    dist = sqrt((x1-y1)^2+(x2-y2)^2);
    
end

function [indices, k] = getNeighbors (data, i, eps )
    % function that returns indices of all neighbors and their number (k)
    k = 0;
    indices = 0;

    for j = 1:size(data,1)
        if i==j
            continue
        end
        if euclidianDist (data(i,1), data(i,2), data(j,1), data(j,2)) <= eps
            k = k+1;
            indices(k) = j;
           
        end
    end
end

function [coreDist] = calcCoreDist (data, neighIndices, i, neighCount, Nmin, coreDistList)
    % function that returns core distance for the selected point (i
    
    % if core distance was already calculated
    if coreDistList(i) ~= -1            % -1 represents undefined
        coreDist = coreDistList(i);
       
    else
    
    % if there is a minimal number of neighbors
        if neighCount >= Nmin-1
        distances = zeros(neighCount,1);
            for k = 1:neighCount
                distances(k) = euclidianDist (data(i,1), data(i,2), data(neighIndices(k),1), data(neighIndices(k),2));
            end   
            distances = sort(distances);
            coreDist = distances(Nmin-1);
        else
        coreDist = -1;
        end
    end
end
   

function [seeds, seedNum, reachDistList] = update(data, neighIndices, neighCount, seeds, seedNum, i, procesList, reachDistList, coreDistList)
%Update the seeds' reachability distance if a smaller value is found.

for k = 1:neighCount
    ind_neigh = neighIndices(k);
    neighbor = data(ind_neigh,:);
    
    if(procesList(ind_neigh)==0) %0 unprocessed
        new_reach = max ( [coreDistList(i); euclidianDist(data(i,1), data(i,2), neighbor(1), neighbor(2))]);
        
        if(reachDistList(ind_neigh)==-1)
            reachDistList(ind_neigh) = new_reach;
            
            seedNum = seedNum +1;
            seeds(seedNum) = ind_neigh;
            
            
        else
            if new_reach < reachDistList(ind_neigh)
                reachDistList(ind_neigh) = new_reach;
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


    
    