clear all; close all; clc

cd("../")
folders=dir("*session*");
if folders(isempty)
    disp("No sessions folder found")
else    
    disp(length(folders)+" session folder found. Now searching for time series volumes data in their respective processed_data folders")
    for f=1:length(folders)
        cd(folders(f).name+"\slug_test\processed_data")
        volumesFile=dir("*volumes.csv");
        if isempty(volumesFile)
            disp("No time series volumes file found. Run imagesToVolumes first to create a time series volumes file from a folder of slug images.")
        else
            T=readtable(volumesFile.name);

            figure
            plot(T.Time, T.Volume_cm3,'xk')
            
            [groupedTimes,groupedVolumes,results]=groupVolumesByTime(T.Time,T.Volume_cm3,100);

        end
    end
end


function [groupedTimes,groupedVolumes,results]=groupVolumesByTime(Times,Volumes,timeGap)
    indices=find((diff(Times)>milliseconds(timeGap))==1);
    for j=1:length(indices)-1
        groupedTimes{j}=Times(indices(j):indices(j+1)-1);     
        groupedVolumes{j}=Volumes(indices(j):indices(j+1)-1);     
        meanVolumes(j)=mean(groupedVolumes{j});
        stdVolumes(j)=std(groupedVolumes{j});
        firstTime(j)=Times(indices(j));
        results=Table(firstTime,meanVolumes,stdVolumes);
    end
end