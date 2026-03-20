clear all; close all; clc

cd("../processed_data")

stdTolerance=50;

volumesFile=dir("*volumes.csv");
if isempty(volumesFile)
    disp("No time series volumes file found. Run imagesToVolumes first to create a time series volumes file from a folder of slug images.")
else
    T=readtable(volumesFile.name);

    figure
    plot(T.time, T.droplet_volume,'xk')

    %Save grouped volumes data and save to "volume_grouped.csv"
    [groupedTimes,groupedVolumes,T2]=groupVolumesByTime(T.time,T.droplet_volume,100);

    %Only consider measurements that have a standard deviation
    %below the threshold
    indices=(T2.std_droplet_volume./T2.mean_droplet_volume*100<stdTolerance);

    T3=table(T2.deposition_time_start(indices),T2.mean_droplet_volume(indices),T2.std_droplet_volume(indices),'VariableNames',["time","mean_droplet_volume","std_droplet_volume"]);
    writetable(T3,"volumes_grouped.csv")

    %Compute volumetric flow rate data and save to "volumetric_flow.csv"
    T3B=table(T2.deposition_time_start(indices),T2.deposition_time_mid(indices),T2.deposition_time_end(indices),T2.mean_droplet_volume(indices),T2.std_droplet_volume(indices),'VariableNames',["deposition_time_start","deposition_time_mid","deposition_time_end","mean_droplet_volume","std_droplet_volume"]);
    volumetricFlowRate=T3.mean_droplet_volume(2:end)/1000./minutes(diff(T3.time));
    corrTimes_start=T3B.deposition_time_start(2:end);
    corrTimes_mid=T3B.deposition_time_mid(2:end);
    corrTimes_end=T3B.deposition_time_end(2:end);

    [V,I]=rmoutliers(volumetricFlowRate);

    T4=table(corrTimes_start(~I),corrTimes_mid(~I),corrTimes_end(~I),V,'VariableNames',["deposition_time_start","deposition_time_mid","deposition_time_end","volumetric_flow_rate"]);
    writetable(T4,"volumetric_flow.csv")
end


function [groupedTimes,groupedVolumes,results]=groupVolumesByTime(Times,Volumes,timeGap)
indices=find((diff(Times)>milliseconds(timeGap))==1);
for j=1:length(indices)-1
    groupedTimes{j}=Times(indices(j):indices(j+1)-1);
    groupedVolumes{j}=Volumes(indices(j):indices(j+1)-1);
    mean_droplet_volume(j)=mean(groupedVolumes{j});
    std_droplet_volume(j)=std(groupedVolumes{j});
    deposition_time_start(j)=groupedTimes{j}(1);
    deposition_time_mid(j)=groupedTimes{j}(round(end/2));
    deposition_time_end(j)=groupedTimes{j}(end);
    results=table(deposition_time_start',deposition_time_mid',deposition_time_end',mean_droplet_volume',std_droplet_volume','VariableNames',["deposition_time_start","deposition_time_mid","deposition_time_end","mean_droplet_volume","std_droplet_volume"]);
end
end