% radi s verzijom 1.4

clear all
instrreset

addpath('C:\Users\Darjan\Desktop\moj_optics\data');
addpath('C:\Users\Darjan\Desktop\moj_optics\functions');
addpath('C:\Users\Darjan\Desktop\moj_optics\functions_help');

s = serial('COM5','BaudRate',115200,'Timeout',500);
s.InputBufferSize = 50000;
fopen(s);

%%
clear data_all

%load('UAE_data_exp_13_complete.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');
%load('UAE_data_exp_09_2018_07_23.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');
load('UAE_data_exp_07_2018_07_16_complete.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');   %dosta sumovit
%load('UAE_data_exp_08_2018_07_19.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');
%load('UAE_data_exp_10_25_07_2018_part_01.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');
%load('UAE_data_exp_12_complete.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');     
%load('UAE_data_exp14_part_03.mat','UAE_single_t_mat_hr','UAE_single_f_mat','UAE_multi_t_mat_hr','UAE_multi_f_mat');

data_all(1,:) = [UAE_single_f_mat; UAE_multi_f_mat(:,1)]/1e3;
data_all(2,:) = [UAE_single_t_mat_hr; UAE_multi_t_mat_hr];

% data_all(1,:) = [UAE_single_f_mat; UAE_multi_f_mat(:,1); UAE_single_f_mat; UAE_multi_f_mat(:,1); UAE_single_f_mat; UAE_multi_f_mat(:,1) ]/2e3;
% data_all(2,:) = [UAE_single_t_mat_hr; UAE_multi_t_mat_hr; UAE_single_t_mat_hr.*-1+2*max(UAE_single_t_mat_hr); UAE_multi_t_mat_hr.*-1+2*max(UAE_multi_t_mat_hr); UAE_single_t_mat_hr+2*max(UAE_multi_t_mat_hr); UAE_multi_t_mat_hr+2*max(UAE_multi_t_mat_hr)];

data_all = sortrows(data_all',2);


%%

Nmin=25;
eps=10;
boxSize=2000;
crosPointNum=500;

large_cluster_perc = 1;
merge_perc = 0.8;
type = 2;
w = 0.5;
t = 160;

firstStage = 0;
iter=0;
compareVect=zeros(size(data_all,1),1);
vrijemeUkupno=0;
vrijemeMin=0;
vrijemeMax=0;

reset=1;

totalStages = floor((size(data_all,1)-crosPointNum)/(boxSize-crosPointNum));
%totalStages =1;
clustMat=zeros(20,100);
errorWindow=zeros(totalStages,1);
reachListStm = zeros(boxSize,1);

figure;
hold on;

pauseTime=0.005;

for step = firstStage:totalStages-1
      
fprintf(s,'%d\n',reset);
pause(pauseTime);

reset=0;

flushinput(s);    
data = data_all(step*(boxSize-crosPointNum)+1:(step+1)*(boxSize)-step*crosPointNum,:);    

pauseTime = 0.005;
for i=1:size(data,1)

    fprintf(s,'%f\n',data(i,1));
    pause(pauseTime);
end

for i=1:size(data,1)

    fprintf(s,'%f\n',data(i,2));
    pause(pauseTime);
end

tic
% optics u matlabu
[orderedList, reachDistList, coreDistList, procesList] = optics(data, Nmin, eps);
 for i = 1:size(data,1)
     orderedReachList(i,1) = reachDistList(orderedList(i));
 end 
[SetClusters, clustnum] = gradient_clustering( orderedReachList, Nmin, t, w, large_cluster_perc, merge_perc, type );
clusterIndices=getClusterIndices(orderedList, SetClusters, clustnum);


clustNum = str2double(fgetl(s)); %
toc
vrijeme=toc;

vrijemeUkupno=vrijeme+vrijemeUkupno;

if vrijeme<vrijemeMin
    vrijemeMin=vrijeme;
end
if step==0
    vrijemeMax=vrijeme;
    vrijemeMin=vrijeme;
end
if vrijeme>vrijemeMax
    vrijemeMax=vrijeme;
end


stmclusterIndices = zeros(boxSize,1);
yMinStart = zeros(clustNum,1);
yMinEnd = zeros(clustNum,1);
yMaxStart = zeros(clustNum,1);
yMaxEnd = zeros(clustNum,1);
timeMin = zeros(clustNum,1);
timeMax = zeros(clustNum,1);

for i=1:size(data,1)
    
    %clusterIndices(i)= fread(s,500,"uint32");
    %clusterIndices(i)= fscanf(s,'%d\n\r');
    stmclusterIndices(i)=str2double(fgetl(s));%
end

for i=1:clustNum
    %clustMat(i,1)=fread(s,clustNum,"uint32");
    %clustMat(i,1)=fscanf(s,'%d\n\r');%
    clustMat(i,step+1)=str2double(fgetl(s));%
end

for i=1:clustNum
   timeMin(i)=str2double(fgetl(s)); 
   timeMax(i)=str2double(fgetl(s)); 
   yMinStart(i)=str2double(fgetl(s));
   yMaxStart(i)=str2double(fgetl(s)); 
   yMinEnd(i)=str2double(fgetl(s)); 
   yMaxEnd(i)=str2double(fgetl(s));   
end

% for i=1:size(data,1)
%     %clusterIndices(i)= fread(s,500,"uint32");
%     %clusterIndices(i)= fscanf(s,'%d\n\r');
%     reachListStm(i)=str2double(fgetl(s));%
% end

colors = ['b','r','g','y','c','m','k'];

for i=1:clustNum 
    fill([timeMin(i) timeMax(i) timeMax(i) timeMin(i) ],[yMinStart(i) yMinEnd(i) yMaxEnd(i) yMaxStart(i) ], colors(mod(clustMat(i,step+1)-1,7)+1));
    ylim([ 50 350]);  
end


% figure
%  for i=1:clustNum 
%         colors = ['b','r','g','y','c','m','k'];
%         index = clusterIndices(:)==i;
%         ylim([ 50 350]);        
%         scatter(data(index,2),data(index,1),'.',colors(mod(clustMat(i,step+1)-1,7)+1)); 
%         %scatter(data(index_ordered,2),data(index_ordered,1),'.');
%         
%         hold on
%  end    
errorWindow(step+1)=numel(find((stmclusterIndices-clusterIndices)~=0));

for i=1:boxSize
    if stmclusterIndices(i)~= clusterIndices(i)
        compareVect(i+step*(boxSize-crosPointNum))=1;
    end
end


end
vrijemeAvg = vrijemeUkupno/totalStages;

% figure 
% scatter(data_all(compareVect==1,2), data_all(compareVect==1,1),'.', 'r');
% hold on
% scatter(data_all(compareVect==0,2), data_all(compareVect==0,1),'.', 'b');

% fclose(s);
% clear s;