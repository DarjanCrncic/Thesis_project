function [timeMin, timeMax, yMinStart, yMinEnd, yMaxStart, yMaxEnd] = feature_extraction(data, clusterIndices, clustNum, timeStampEnd, timeStampStart, crosPointNum, newClustTest, prevClustVect, currClustVect, prevyMaxEnd, prevyMinEnd, Nmin)
N = size(data, 1);

timeMin = zeros(clustNum,1);
timeMax = zeros(clustNum,1);

yMinStart = zeros(clustNum,1);
yMinEnd = zeros(clustNum,1);
yMaxStart = zeros(clustNum,1);
yMaxEnd = zeros(clustNum,1);
minTpoint = zeros(clustNum,1);
maxTpoint = zeros(clustNum,1);
 
for i=1:N
     index = clusterIndices(i);
     
     if index ~= 0
         if timeMin(index) == 0
             timeMin(index) = data(i,2);
             minTpoint(index) = i;
         end
         timeMax(index)=data(i,2);
         maxTpoint(index) = i;  
         
         if data(i,1)> yMaxEnd(index)
             yMaxEnd(index)=data(i,1);
         end
         if data(i,1)< yMinEnd(index) || yMinEnd(index) == 0
             yMinEnd(index)=data(i,1); 
         end  
     end
end

 for i=1:clustNum    
     
    % cluster started in this window
    if newClustTest(i) == 1 
        
        yMinStart(i) = 0;
        yMaxStart(i) = 0;
        counter = 0;
        start = 1;
        
        while(counter<Nmin && start<=N )
            index = clusterIndices(start);
            if(index == i)
                counter=counter+1;
                if data(start,1) < yMinStart(i) || yMinStart(i)==0 
                    yMinStart(i)=data(start,1);
                end
                if data(start,1) > yMaxStart(i) 
                    yMaxStart(i) = data(start,1);
                end               
            end
            start=start+1;
        end
    end
    
    % if cluster ends in this window
    if timeMax(i)-timeStampStart < 0.81*(timeStampEnd-timeStampStart)
        
        yMinEnd(i) = 0;
        yMaxEnd(i) = 0;
        counter = 0;
        start = N;

        while(counter<Nmin && start>0)
            index = clusterIndices(start);
            if(index == i)
                counter=counter+1;
                if data(start,1) < yMinEnd(i) || yMinEnd(i)==0 
                    yMinEnd(i)=data(start,1);
                end
                if data(start,1) > yMaxEnd(i) 
                    yMaxEnd(i) = data(start,1);
                end               
            end
            start=start-1;
        end
    else
        timeMax(i) = timeStampEnd;
    end
 
    
    % if its an older cluster, search for end points in previous iteration
    if newClustTest(i) == 0
       timeMin(i)=timeStampStart;
       
       if timeMax(i) < timeStampStart
          timeMax(i) = timeStampStart;          
       end
          
       index = currClustVect(i);
       j=1;
       
       while(prevClustVect(j)~=0)
           if prevClustVect(j)==index
              yMinStart(i)=prevyMinEnd(j);
              yMaxStart(i)=prevyMaxEnd(j);
              break;
           end
           j=j+1;
       end
    end   
 end
 

 
 
end
