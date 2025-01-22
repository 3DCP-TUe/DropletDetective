clear all; close all; clc;

%Optionally give the main folder to analyse. If not specified, the folder
%in which the script is located will be used.
cd("D:\Chapter 9 Slug images")

% Give the range of the region of interest
ROIylim=[307 483];
ROIxlim=[180 1329];

% Give the minimum size of a droplet in pixels
minSizeXDir=100;
minSizeYDir=10;

%% Calibration
% cd("Calibration")
% file=dir("*.jpg");
%
% imshow(imread(file.name));
% cd("../")

% The characteristic length is determined manually in pixels
length_pix=470-149; %pix
length_mm=119.67; %mm
calibration=length_mm/length_pix; %%mm/pix

%%
folders=dir("*slug_test_images");

for f=1:length(folders)
    folder=folders(f).name;
    cd(folder)

    files=dir("*.jpg");
    
    %Create empty arrays
    dur=duration(strings(length(files),1),'Format','hh:mm:ss.SSS');
    volume=zeros(length(files),1);

    %Loop over all of the images
    for i=1:length(files)
        image=imread(files(i).name);
        ROI=image(ROIylim(1):ROIylim(2),ROIxlim(1):ROIxlim(2),:);
        
        %Show the region of interest once, to allow for a check
        if i==1
            figure
            imshow(ROI)
        end

        %Only consider dark regions
        filteredROI=ROI(:,:,1)>20;

        %Detect objects by checking every column of pixels to see if it
        %contains more than minSizeYDir black pixels
        detectVert=sum(filteredROI,1)<diff(ROIylim)-minSizeYDir;

        indicesA=find(diff(detectVert(1:end))==1);
        indicesB=find(diff(detectVert(1:end))==-1);

        if length(indicesA)>1
            %Only analyse the largest object
            [~,I]=max(diff([indicesA(1:length(indicesB)); indicesB],1));

            filteredROI2=abs(filteredROI-1);
            volume(i)=sum((sum(filteredROI2(:,indicesA(I):indicesB(I)),1)./2).^2*pi)*calibration^3; %mm³
        else
            volume (i)=NaN;
        end
        
        %Show the filtered ROI once, to allow for a check
        if i==1
            figure
            imshow(filteredROI2(:,indicesA(I):indicesB(I)))
        end

        %Obtain the timestamp from the image title
        timetemp=str2double(split(convertCharsToStrings(files(i).name(14:end-4)),"."));
        dur(i)=duration(timetemp(1),timetemp(2),timetemp(3),timetemp(4),'Format','hh:mm:ss.SSS');
        if round(i/1000)==i/1000
            disp(i/length(files)*100+"%")
        end
    end
    cd("../")

    %Write the results to a timeseries csv file
    T=table(dur,volume/1000,'VariableNames',["Time","Volume_cm3"]);
    writetable(T,folder+"_processed.csv");
end