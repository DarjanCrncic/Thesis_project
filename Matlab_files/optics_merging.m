function [clustersMatrix, prevCluster, newClustInd, newClustTest, colorMatrix] = optics_merging(clusterIndices, clustersMatrix, prevCluster,  crosPointNum, clustnum, newClustInd, iter, step, colorMatrix )

newClustTest = zeros(clustnum,1);

if iter>0   
    currCluster = zeros(crosPointNum,1);
    simMatrix = zeros(100,100);
 
% labeling first "crosPointNum" amount of points depending on their cluster 
    for i=1:crosPointNum    
        currCluster(i) = clusterIndices(i);
    end
         
%filling simmilarity matrix
    for i=1:crosPointNum
        if prevCluster(i)~=0 && currCluster(i)~=0
            simMatrix(prevCluster(i), currCluster(i))=simMatrix(prevCluster(i), currCluster(i))+1;    
        end
    end 
     

% combining clusters
    if clustnum>0
        for i=1:clustnum  
                pom = simMatrix(:,i);              
                [value, index] = max(pom);
                pom2 = simMatrix(index,:);
                
%                 if value >= 0.8*sum(pom) && value >= 0.8*sum(pom2) && value~=0
%                     clustersMatrix(i,step+1)=clustersMatrix(index,step);
%                     colorMatrix(i,step+1)=colorMatrix(index,step);
%                 elseif value >= 0.8*sum(pom) && value~=0
%                     clustersMatrix(i, step+1)=newClustInd;
%                     colorMatrix(i,step+1)=clustersMatrix(index,step);
%                     newClustInd=newClustInd+1;
%                     newClustTest(i)=1;
%                 else
%                     clustersMatrix(i, step+1)=newClustInd;
%                     if value==0
%                         colorMatrix(i,step+1)=newClustInd;
%                     else
%                         colorMatrix(i,step+1)=colorMatrix(index,step);
%                     end
%                     newClustInd=newClustInd+1;
%                     newClustTest(i)=1;
%                     
%                 end
               
                if value >= 0.7*sum(pom) && value >= 0.7*sum(pom2) && value~=0
                    clustersMatrix(i,step+1)=clustersMatrix(index,step);
                else               
                    clustersMatrix(i, step+1)=newClustInd;
                    newClustInd=newClustInd+1;
                    newClustTest(i)=1;
                end
        end
    end
end

if iter==0
    for i=1:clustnum
        clustersMatrix(i, 1)=i;
        colorMatrix(i,1)=i;
        newClustTest(i)=1;
    end
    newClustInd = clustnum+1;
    

end

% saving "crosPointNum" amount of cluster labels

 prevCluster = zeros(crosPointNum,1);

 for i=1:crosPointNum     
     prevCluster(i)= clusterIndices(i+length(clusterIndices)-crosPointNum);     
 end

end