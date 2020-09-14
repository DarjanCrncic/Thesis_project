
function [SetClusters, clustnum] = gradient_clustering( reachDistList, MinPoints, t, w, bound_perc, merge_perc, minMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for gradient clustering of resulting data from optics algorithm

%INPUTS
% - reachDistList (reachability distance list, MUST BE CORRECTLY INDEXED)
% - MinPoints (minimal number of points to form a cluster, same as Nmin)
% - t (angle of steepness for checking inflections, recomended 160 )
% - w (distance between points, recomended 0.5)
% - bound_perc (maximum cluster size, number from 0 to 1)
% - merge_perc (similarity factor for merging clusters)
% - minMax (minMax==1 keeps only the smallest clusters, minMax==2 keeps the
%           biggest ones)

%OUTPUTS
% - SetClusters (final result, each row presents a cluster, first number is     
%                the cluster start point, second number cluster end point, 
%                third number cluster size)
% - clustNum (number of clusters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(reachDistList)
 if reachDistList(i) == -1 
     reachDistList(i) = 1024;  
 end
end

t = cos(t*pi/180);
ListSize = size(reachDistList,1);
LastEnd = ListSize;

CurrCluster = 0;
clustnum = 0;

stpnum = 1;
StartPoints = zeros(LastEnd,1);
StartPoints(stpnum) = 1;

for i = 2:(size(reachDistList,1)-1)                                     
    if inflectionIndex(reachDistList, i, w) > t                             %ako je tocka infleksija
        if gradientDeterminant(reachDistList, i, w) > 0                     %ako je gradijent pozitivan
            
            if length(CurrCluster) >= MinPoints                             %ako je trenutni cluster dovoljno velik, spremi ga
                clustnum = clustnum+1;
                SetOfClusters{clustnum,1} = CurrCluster;
            end
            
            CurrCluster = 0;
            
            if ~isempty(StartPoints)                                        %ako ima startnih tocaka
                
                if reachDistList(StartPoints(stpnum)) <= reachDistList(i)
                    StartPoints=StartPoints(1:stpnum-1);
                    stpnum = stpnum-1;
                end
            end
            
            if ~isempty(StartPoints)
                
                while(reachDistList(StartPoints(stpnum)) < reachDistList(i))
                    
                    TempCluster = StartPoints(stpnum):LastEnd;
                    
                    if length(TempCluster) >= MinPoints
                        clustnum = clustnum+1;
                        SetOfClusters{clustnum,1} = TempCluster;
                    end
                    
                    stpnum = stpnum-1;
                    StartPoints=StartPoints(1:stpnum);
                   
                    
                 end
                
                TempCluster = StartPoints(stpnum):LastEnd;
                
                if length(TempCluster) >= MinPoints
                        clustnum = clustnum+1;
                        SetOfClusters{clustnum,1} = TempCluster;
                end
            end
            
            if reachDistList(i+1) < reachDistList(i)
                stpnum = stpnum+1;
                StartPoints(stpnum) = i;
            end
            
        else
                
            if reachDistList(i+1) > reachDistList(i) %&&  reachDistList(i)>=reachDistList(stpnum)
                  LastEnd = i+1;
                   if stpnum == 0
                       CurrCluster = 1:LastEnd;
                  
                   else
                  CurrCluster = StartPoints(stpnum):LastEnd;
                  end
            end
                
        end
    end
end



while (~isempty(StartPoints))
                
     CurrCluster = StartPoints(stpnum): size(reachDistList,1);
        
     if (reachDistList(StartPoints(stpnum)) > reachDistList(size(reachDistList,1))) && (length(CurrCluster) >= MinPoints)
            clustnum = clustnum+1;
            SetOfClusters{clustnum,1} = CurrCluster;
     end
     
     StartPoints=StartPoints(1:stpnum-1);
     stpnum = stpnum-1;
end


%%
if clustnum==0
    SetClusters = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clustnum>=1 %%POSTPROCESSING!!!!!!!!!!!!!!
    
for i=1:clustnum
   setClustersStart(i)=SetOfClusters{i,1}(1); 
   setClustersEnd(i)=SetOfClusters{i,1}(end); 
end
    
   
%% 1. removing clusters detected cos of sensor resolution%%%%%%%%%%

%    i=1;
%    while( i<=clustnum )
%        if reachDistList(setClustersStart(i))<=1  || reachDistList(setClustersEnd(i)) <=1 %depends on scaling number!!!
%         clustnum = clustnum-1;
%         setClustersStart(i)=[];
%         setClustersEnd(i)=[];
%         continue;
%        end
%       
%        i=i+1;
%        
%    end
%%%%%%%%%%%%%%%%%%%%%

% figure
% plot(reachDistList);
% hold on
% axis([0 size(reachDistList,1) 0 max(reachDistList(reachDistList~=1024))*3]);
% set(gca,'Box','on');
%  title('Dohvatna krivlulja, set podataka 8');
%  ylabel('dohvatna udaljenost');
%  xlabel('redoslijed obraðivanih toèaka');
% 
%  for i=1:clustnum
% 
%     first = setClustersStart(i);
%     last = setClustersEnd(i);
%     line([first;last],[reachDistList(first);reachDistList(last)],'color',rand(1,3));
% 
%  end 

%% 2. geting rid of large clusters, larger than bound_perc points

    i = 1;
    while ( i<=clustnum ) 
        size_clust = setClustersEnd(i)-setClustersStart(i);
        if size_clust >= bound_perc*size(reachDistList,1) 
            clustnum = clustnum-1;
            setClustersStart(i)=[];
            setClustersEnd(i)=[];
            continue;
        end
        if reachDistList(setClustersStart(i)) == 1024
            setClustersStart(i)=setClustersStart(i)+1;    
            continue;
        end
        if reachDistList(setClustersEnd(i)) == 1024           
            setClustersEnd(i)=setClustersEnd(i)-1; 
            continue;
        end
            i = i+1;                  
    end    
    
% figure
% plot(reachDistList);
% hold on
% axis([0 size(reachDistList,1) 0 max(reachDistList(reachDistList~=1024))*3]);
% set(gca,'Box','on');
% title('Dohvatna krivlulja, prikaz otkrivenih grupacija');
% ylabel('dohvatna udaljenost');
% xlabel('redoslijed obraðivanih toèaka');
%  grid on
%  ylim([0 20]);
%  for i=1:clustnum
% 
%     first = setClustersStart(i);
%     last = setClustersEnd(i);
%     line([first;last],[reachDistList(first);reachDistList(last)],'color',rand(1,3));
% 
%  end 


%% 3. removing clusters that spread over two clusters divided by noise point  
i=1;
while(i<=clustnum)
       wrong=0;       
      
       for j=setClustersStart(i):setClustersEnd(i)
          if reachDistList(j)==1024
              wrong=1;
              break;
          end
       end
         
       if wrong==1
           setClustersStart(i)=[];
           setClustersEnd(i)=[];
           clustnum=clustnum-1;
           continue
       end
       i=i+1;
end 

% figure
% plot(reachDistList);
% hold on
% axis([0 size(reachDistList,1) 0 max(reachDistList(reachDistList~=1024))*3]);
% set(gca,'Box','on');
% title('Removing clusters 3. step');
% ylabel('Reachability distance');
% xlabel('Ordered index');
% 
%  for i=1:clustnum
% 
%     first = setClustersStart(i);
%     last = setClustersEnd(i);
%     line([first;last],[reachDistList(first);reachDistList(last)],'color',rand(1,3));
% 
%  end 
 

%% 4. merging simmilar clusters
 
% sorting clusters by size
    for i=2:clustnum
            value = setClustersEnd(i)-setClustersStart(i);
            mem1 = setClustersStart(i);
            mem2 = setClustersEnd(i);
            hole = i;

            while ( hole>1 && setClustersEnd(hole-1)-setClustersStart(hole-1)<value )
                setClustersStart(hole) =  setClustersStart(hole-1);
                setClustersEnd(hole) =  setClustersEnd(hole-1);
                hole=hole-1;
            end
            setClustersStart(hole) =  mem1;
            setClustersEnd(hole) =  mem2;
    end

% merging
 i=1;
 while(i<clustnum)
            j=i+1;
            while(j<=clustnum)
                if setClustersEnd(i)>=setClustersEnd(j) && setClustersStart(i)<=setClustersStart(j) 
                    clustSize=setClustersEnd(i)-setClustersStart(i);

                    if setClustersEnd(j)-setClustersStart(j)>=clustSize*merge_perc									
                        for k=j:clustnum-1
                            setClustersStart(k)=setClustersStart(k+1);
                            setClustersEnd(k)=setClustersEnd(k+1);
                        end 
                            setClustersStart(clustnum)=[];
                            setClustersEnd(clustnum)=[];
                            clustnum=clustnum-1;
                            continue;  

                    end                    
                end
                j=j+1;
            end
            i=i+1;
 end

% figure
% plot(reachDistList);
% hold on
% axis([0 size(reachDistList,1) 0 max(reachDistList(reachDistList~=1024))*3]);
% set(gca,'Box','on');
% title('Removing clusters 4. step');
% ylabel('Reachability distance');
% xlabel('Ordered index');
% 
%  for i=1:clustnum
% 
%     first = setClustersStart(i);
%     last = setClustersEnd(i);
%     line([first;last],[reachDistList(first);reachDistList(last)],'color',rand(1,3));
% 
%  end 


%% 5. minMax==1 keeps only the smallest clusters, minMax==2 keeps the biggest ones
if minMax==1
end
%%% možda ni neæu implementirat

if minMax==2

%keeping only the largest clusters
    i=1;
    while(i<clustnum)
        j=i+1;
        while(j<=clustnum)        
            if setClustersStart(j)>=setClustersStart(i) && setClustersEnd(j)<=setClustersEnd(i)
                setClustersStart(j)=[];
                setClustersEnd(j)=[];
                clustnum=clustnum-1;
                continue;
            end
            j=j+1;
        end
        i=i+1;
    end
    
end 


%%
for i=1:clustnum
    setClustersEnd(i)=setClustersEnd(i)-1;
end

if clustnum>0
    SetClusters=zeros(clustnum,3);
    for i=1:clustnum
        SetClusters(i,1)=setClustersStart(i);
        SetClusters(i,2)=setClustersEnd(i);
        SetClusters(i,3)=SetClusters(i,2)-SetClusters(i,1); 
    end
else
    SetClusters=0;
end

% colors = ['b','r','g','y','c','m','k'];
% figure
% plot(reachDistList,'linewidth',1);
% hold on
% axis([0 size(reachDistList,1) 0 max(reachDistList(reachDistList~=1024))*3]);
% set(gca,'Box','on');
% ylim([0 20]);
% grid on
% title('Rezultat gradijentnog grupiranja');
% ylabel('dohvatna udaljenost');
% xlabel('redoslijed obraðivanih toèaka');
% set(gca,'FontSize',14)
% 
%  for i=1:clustnum
% 
%     first = setClustersStart(i);
%     last = setClustersEnd(i);
%     line([first;last],[reachDistList(first);reachDistList(last)],'color',rand(1,3));
% 
%  end 



end
end      

function [vectAbs] = reachabilityVector (rx, ry, w)

vectAbs = sqrt((ry-rx)^2 + w^2);


end

function [inflecInd] = inflectionIndex( reachDistList, y, w)
 % calculate inflection index for point y 
 % w is distance between points in reachDistList (1?)
 
 x_r = reachDistList(y-1);
 y_r = reachDistList(y);
 z_r = reachDistList(y+1);
 
 prev_vector = reachabilityVector(x_r,y_r,w);
 next_vector = reachabilityVector(y_r,z_r,w);
 
 inflecInd = (-w^2 + (x_r-y_r)*(z_r-y_r))/(prev_vector*next_vector);
 
end

function [gradDet] = gradientDeterminant( reachDistList, y, w)

 x_r = reachDistList(y-1);
 y_r = reachDistList(y);
 z_r = reachDistList(y+1);

 gradDet = w*(y_r-x_r) - w*(z_r-y_r);
 
end



