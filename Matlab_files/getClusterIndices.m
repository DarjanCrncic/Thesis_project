function [clusterIndices] = getClusterIndices(orderedList, SetClusters, clustNum)
    
    listSize = length(orderedList);
    clusterIndices = zeros(listSize,1);
    clusterIndicesReal = zeros(listSize,1);
    
    for i=1:listSize
        for j=1:clustNum
            index = orderedList(i);
            if i >= SetClusters(j,1) && i <= SetClusters(j,2)
                clusterIndices(index) = j;
               
            end
        end
    end
      

end