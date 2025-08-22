%{
This file is based on Droplet Detective. Droplet Detective is licensed under the terms 
of GNU General Public License as published by the Free Software Foundation. For more 
information and the LICENSE file, see <https://github.com/3DCP-TUe/DropletDetective>.

This file should be placed in a folder called "scripts". It automatically
searches one folder level up for folders containing the keyword "session".
In those folders, it opens a folder called "slug_test" and in there, it
opens the raw_data folder. It processes this data depending on the type of
raw data present.
%}


clear all; close all; clc;

%Navigate to main folder
cd("../")

%Search for sessions
folders=dir("*session*");
if isempty(folders)
    disp("No sessions folder found")
else    
    disp(length(folders)+" session folders found.")
    for f=1:length(folders)
        cd(folders(f).name)
        disp("--------------------")
        disp("Opening: "+folders(f).name)
        disp("--------------------")
        if not(isfolder("slug_test"))
            disp("WARNING: No folder named slug_test is found.")
        else
            cd("slug_test")
            if not(isfolder("raw_data"))
                disp("WARNING: No folder named raw_data is found.")
            else
                cd("raw_data")
                files=dir("*.csv");
                cd("../")
            end

            if not(isfolder("processed_data"))
                mkdir("processed_data")
            end
                cd("processed_data")

                %Analyze load cell
                if exist(files(1).name(1:end-4)+"_processed_yield_stress.csv")~=2 && exist(files(1).name(1:end-4)+"_processed_mass_flow.csv")~=2
                    analysisLoadCell("py",25,12,true)
                else
                    disp("Load cell data not analyzed since the following files already exist in the processed_data folder:")
                    disp(files(1).name(1:end-4)+"_processed_yield_stress.csv")
                    disp(files(1).name(1:end-4)+"_processed_mass_flow.csv")
                end
                
                %Analyze images
                if exist(files(1).name(1:end-4)+"_processed_volumes.csv")~=2
                    % Calibration
                    % cd("Calibration")
                    % file=dir("*.jpg");
                    %
                    % CAL=imread(file.name);
                    % filteredCAL=CAL(:,:,1)>20;
                    % imshow(filteredCAL)
                    % cd("../")

                    % The characteristic length is determined manually in pixels
                    length_pix=470-149; %pix
                    length_mm=119.67; %mm
                    calibration=length_mm/length_pix; %%mm/pix
                    imagesToVolumes([180 1329],[307 483],10,0.3728)
                else
                    disp("Slug image data not analyzed since the following files already exist in the processed_data folder:")
                    disp(files(1).name(1:end-4)+"_processed_volumes.csv")
                end

                %Analyze volumes
                if exist(files(1).name(1:end-4)+"_processed_volumes_grouped.csv")~=2 && exist(files(1).name(1:end-4)+"_processed_volumetric_flow.csv")~=2
                    analyzeVolumes(100,50)
                else
                    disp("Slug image data not analyzed since the following files already exist in the processed_data folder:")
                    disp(files(1).name(1:end-4)+"_processed_volumes_grouped.csv")
                    disp(files(1).name(1:end-4)+"_processed_volumetric_flow.csv")
                end


                cd("../")
                cd("../")
                disp("-----FINISHED-----")
        end
        cd("../")
    end


end


function analysisLoadCell(loggerType,nozzleDiameter,bucketWeight,plotting)

%% File locations and settings

% Select the loggerType
% "py" refers to the OPC-UA logger in python https://github.com/arjendeetman/Python-OPC-UA
% "ua" refers to the logger in UA Expert

% -------------------------------------------------------------------------
% Set values for stable plateau detection
% -------------------------------------------------------------------------

% A moving standard deviation is applied with window size "k1". When the
% value of this moving standard deviation is larger than "lim1" that
% measured value is labelled as unstable.

% Window size for movind standard deviation
k1 = 10;

% Tolerance bound
lim1 = 0.05; %N

% When the time inbetween two stable measurments is larger than "interTime",
% the stable measurements belong to a different stable plateau.
interTimeLimit = 0.2; %seconds

% Minimum number of stable measurements to select the stable plateau for
% further processing
minNo = 5;

% -------------------------------------------------------------------------
% Select the method that is used for filtering. Options are:

% "removeOutliers": uses the standard rmoutliers method from Matlab

% "medianIntervalTime": removes outliers based on their distance from the
% moving median of the time between droplets.

% "medianMass": removes outliers based on their distance from the
% moving median of the mass of the droplets.

% "medianMass+stdMass": removes outliers based on their distance from the
% moving median of the mass of the droplets and subsequently removes
% outliers based on the an "factorStd" number of (moving) standard
% deviations away from the median.
% -------------------------------------------------------------------------

filterMethod = "medianMass+stdMass";

% -------------------------------------------------------------------------
% Select tolerance and movingInteger for filterMethod "medianIntervalTime",
% "medianMass", and "medianMass+stdMass".

% tolerance: Upper (1+tolerance) and lower (1-tolerance) limit of the allowed tolerance
% as a factor of the % median interval time (for "medianIntervalTime") or
% mass (for "medianMass", and "medianMass+stdMass").
% Recommended value between 0.5 and 1.

% movingInteger: Window size over which the moving median or moving standard
% deviation is calculated.
% ------------------------------------------------------------------------

tolerance=0.8;
movingInteger = 200;

% -------------------------------------------------------------------------
% Settings for filterMethod "medianMass+stdMass"

% Tolerance bound defined by the number of standard deviations that are
% allowed (for example: 2 sigma).
% -------------------------------------------------------------------------

factorStd=2;

%% No more settings below this line.
%% Read file
T=[];
if exist('filePath')==1
    cd(filePath)
    T=[T; readtable(fileName+fileExtension,'Delimiter',',')];
else
    cd("../raw_data")
    files=dir("*.csv");
    disp("Analyzing file: "+files(1).name)
    for i=1:length(files)
        T=[T; readtable(files(i).name,'Delimiter',',')];
    end
end


if loggerType=="ua"
    for j=1:height(T)
        Ttemp=datetime(T.SourceTimeStamp{j},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','UTC');
        T2(j,1)=duration(hour(Ttemp)-12,minute(Ttemp),second(Ttemp));
    end
elseif loggerType=="py"
    T2=T.Time;
    T.Value=T.Load;
end

%% Detect slugs
Dur=T2-T2(1);

Val1=T.Value(T.Value>bucketWeight);
Dur1=Dur(T.Value>bucketWeight);
Time1=T2(T.Value>bucketWeight);

S=movstd(Val1,k1);

Val2=Val1(S<lim1); %Stable measurements
Dur2=Dur1(S<lim1); %Duration corresponding to stable measurements
Time2=Time1(S<lim1); %Times corresponding to stable measurements

clear S Val3 Dur3 firstDur lastDur meanLoad firstTime lastTime
%This loops over all stable measurements
k=1;
j=0;
for i=2:length(Val2)
    %   diffDur(i-1)=seconds(Dur2(i)-Dur2(i-1));
    j=j+1;
    if Dur2(i)-Dur2(i-1)>seconds(interTimeLimit)
        if j<minNo
            Val3{k}=0;
            Dur3{k}=NaT-NaT;
            j=1;
        else
            meanLoad(k)=mean(Val3{k});
            firstDur(k)=Dur3{k}(1);
            lastDur(k)=Dur3{k}(end);
            firstTime(k)=Time3{k}(1);
            lastTime(k)=Time3{k}(end);

            k=k+1;
            j=1;
        end
    end
    Val3{k}(j)=Val2(i);
    Dur3{k}(j)=Dur2(i);
    Time3{k}(j)=Time2(i);
end

% Calculate mass flow per bucket
clear meanLoad2 firstTime2
k=1;
j=1;
for i=2:length(meanLoad)
    if meanLoad(i)>meanLoad(i-1)-1 && diff([firstTime(i-1) firstTime(i)])<seconds(10)
        meanLoad2{k}(j)=meanLoad(i-1);
        firstTime2{k}(j)=firstTime(i-1);
        j=j+1;
    else
        j=1;
        k=k+1;
    end
end

figure
raw=plot(T2,T.Value);
hold on
for k2=1:length(meanLoad2)
    plot(firstTime2{k2},meanLoad2{k2})
end

k1=0;
for k2=1:length(meanLoad2)
    if isempty(meanLoad2{k2})==0
        k1=k1+1;
        % massFlowBucket(k1)=(meanLoad2{k2}(end)-meanLoad2{k2}(1))/9.81/minutes(firstTime2{k}(end)-firstTime2{k}(1));
        midTimeBucket(k1)=mean([firstTime2{k2}(end); firstTime2{k2}(1)]);
        startTimeBucket(k1)=firstTime2{k2}(1);
        endTimeBucket(k1)=firstTime2{k2}(end);
        x=minutes(firstTime2{k2}'-firstTime2{k2}(1));
        b1=([ones(length(x),1) x])\(meanLoad2{k2}'/9.81);
        massFlowBucket(k1)=b1(2);
    end
end
%
% k1=1;
% k2=1;
% figure
% x=minutes(firstTime2{k2}'-firstTime2{k2}(1));
% y=(meanLoad2{k2}'/9.81);
% yCalc1 = massFlowBucket{k1}(1)+massFlowBucket{k1}(2)*x;
% scatter(x,y)
% hold on
% plot(x,yCalc1)

massFlowA=diff(meanLoad/9.81)./minutes(diff(firstDur));
massFlow=massFlowA(massFlowA>0.1*median(massFlowA)&massFlowA<1.9*median(massFlowA)); %kg/min
%%
figure
plot(midTimeBucket,massFlowBucket)
%%

%Filter slug values based on filterMethod.
if filterMethod =="removeOutliers"
    [SlugMass, I]=rmoutliers(diff(meanLoad)');
    firstDur2a=firstDur(2:end);
    firstTime2a=firstTime(2:end);
    ContactTimeCorr=firstDur2a(~I)';
    ContactTime=firstTime2a(~I)';
    disp("Outliers removed: "+sum(I)+" of "+length(SlugMass))
elseif filterMethod == "medianIntervalTime"
    interTime=seconds(firstTime(2:end)-lastTime(1:end-1));
    SlugMassA=diff(meanLoad)';
    SlugMassB=SlugMassA(SlugMassA>0.1*median(SlugMassA));

    interTimeB=interTime(SlugMassA>0.1*median(SlugMassA));
    movMedianInterTimeB=movmedian(interTimeB,movingInteger);

    SlugMass=SlugMassB(interTimeB>(1-tolerance)*movMedianInterTimeB&interTimeB<(1+tolerance)*movMedianInterTimeB);

    firstDur2a=firstDur(2:end);
    ContactTimeCorrA=firstDur2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTimeCorr=ContactTimeCorrA(interTimeB>(1-tolerance)*movMedianInterTimeB&interTimeB<(1+tolerance)*movMedianInterTimeB);

    firstTime2a=firstTime(2:end);
    ContactTimeA=firstTime2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTime=ContactTimeA(interTimeB>(1-tolerance)*movMedianInterTimeB&interTimeB<(1+tolerance)*movMedianInterTimeB);
elseif filterMethod == "medianMass"
    SlugMassA=diff(meanLoad)';
    SlugMassB=SlugMassA(SlugMassA>0.1*median(SlugMassA));

    movMedianMass=movmedian(SlugMassB,movingInteger);

    SlugMass=SlugMassB(SlugMassB>(1-tolerance)*movMedianMass&SlugMassB<(1+tolerance)*movMedianMass);

    firstDur2a=firstDur(2:end);
    ContactTimeCorrA=firstDur2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTimeCorr=ContactTimeCorrA(SlugMassB>(1-tolerance)*movMedianMass&SlugMassB<(1+tolerance)*movMedianMass);

    firstTime2a=firstTime(2:end);
    ContactTimeA=firstTime2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTime=ContactTimeA(SlugMassB>(1-tolerance)*movMedianMass&SlugMassB<(1+tolerance)*movMedianMass);
elseif filterMethod == "medianMass+stdMass"
    SlugMassA=diff(meanLoad)';
    SlugMassB=SlugMassA(SlugMassA>0.1*median(SlugMassA));

    SlugMassC=SlugMassB(SlugMassB>(1-tolerance)*movmedian(SlugMassB,movingInteger)&SlugMassB<(1+tolerance)*movmedian(SlugMassB,movingInteger));
    SlugMass=SlugMassC(SlugMassC>movmedian(SlugMassC,movingInteger)-factorStd*movstd(SlugMassC,movingInteger)&SlugMassC<movmedian(SlugMassC,movingInteger)+factorStd*movstd(SlugMassC,movingInteger));

    firstDur2a=firstDur(2:end);
    ContactTimeCorrA=firstDur2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTimeCorrB=ContactTimeCorrA(SlugMassB>(1-tolerance)*movmedian(SlugMassB,movingInteger)&SlugMassB<(1+tolerance)*movmedian(SlugMassB,movingInteger));
    ContactTimeCorr=ContactTimeCorrB(SlugMassC>movmedian(SlugMassC,movingInteger)-factorStd*movstd(SlugMassC,movingInteger)&SlugMassC<movmedian(SlugMassC,movingInteger)+factorStd*movstd(SlugMassC,movingInteger));

    firstTime2a=firstTime(2:end);
    ContactTimeA=firstTime2a(SlugMassA>0.1*median(SlugMassA))';
    ContactTimeB=ContactTimeA(SlugMassB>(1-tolerance)*movmedian(SlugMassB,movingInteger)&SlugMassB<(1+tolerance)*movmedian(SlugMassB,movingInteger));
    ContactTime=ContactTimeB(SlugMassC>movmedian(SlugMassC,movingInteger)-factorStd*movstd(SlugMassC,movingInteger)&SlugMassC<movmedian(SlugMassC,movingInteger)+factorStd*movstd(SlugMassC,movingInteger));
end

YieldStress=SlugMass/(sqrt(3)*1/4*pi*nozzleDiameter^2)*1000; %kPa

%% Display mass flow rate

disp("Sum "+sum(SlugMass)/9.81+" kg - "+length(SlugMass) +" slugs")
disp(mean(massFlow)+" kg/min")


%% Plot and save results
cd("../")
if not(isfolder("processed_data"))
    mkdir("processed_data")
end
cd("processed_data")
if plotting == true

    if not(isfolder("figures"))
        mkdir("figures")
    end
    cd("figures")

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 500 1920 500];
    raw=plot(Dur+T2(1),T.Value);
    hold on
    % xlim([minutes(1),minutes(2)])
    stable=plot(Dur2+T2(1),Val2,'.r');
    % firstDurPlot=plot(firstDur,meanSlug,'xg');
    % lastDurPlot=plot(lastDur,meanSlug,'xb');
    yyaxis right
    hold on
    pl=plot(ContactTime,SlugMass,'x','Color','m');
    ax=gca;
    ax.YAxis(2).Color='m';

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 250 1920 500];
    grid on
    box on
    plot(firstDur2a+T2(1),diff(meanLoad),'xk')
    hold on
    plot(ContactTime,SlugMass,'xr')
    ylabel('Mass [N]')

    fig=figure;
    fig.Units='centimeters';
    fig.Position=[10 10 30 9];
    hold on
    grid on
    box on
    plot(ContactTime,YieldStress,'xk')
    xlabel('Time','Interpreter','latex','FontSize',14)
    ylabel('Yield stress [kPa]','Interpreter','latex','FontSize',14)
    saveas(fig,files(1).name(1:end-4)+"_yield_stress_time.svg")

    fig=figure;
    fig.Units='centimeters';
    fig.Position=[10 10 14 9];
    hold on
    grid on
    box on
    histogram(YieldStress,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5)
    xlabel('Yield stress [kPa]','Interpreter','latex','FontSize',14)
    ylabel('Frequency','Interpreter','latex','FontSize',14)
    saveas(fig,files(1).name(1:end-4)+"_processed_histogram.svg")
    cd("../")
end

%% Save as csv file
saveData=table(ContactTime,SlugMass,YieldStress,'VariableNames',["time","droplet_mass","yield_stress"]);
writetable(saveData,files(1).name(1:end-4)+"_processed_yield_stress.csv")

saveData2=table(startTimeBucket',endTimeBucket',midTimeBucket',massFlowBucket','VariableNames',["collection_time_start","collection_time_mid","collection_time_end","mass_flow"]);
writetable(saveData2,files(1).name(1:end-4)+"_processed_mass_flow.csv")
end

function imagesToVolumes(ROIxlim,ROIylim,minSizeYDir,calibration)

cd("../raw_data")

%%
folders=dir("*slug_test_images");
if isempty(folders)
    disp("No folders found. Folders ending in slug_test_images should be unzipped before running this script")
else
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
        cd("../../processed_data")

        %Write the results to a timeseries csv file
        T=table(dur,volume/1000,'VariableNames',["time","droplet_volume"]);
        writetable(T,folder(1:8)+"_inline_"+f+"_processed_volumes.csv");
        
    end
end

end

function analyzeVolumes(timeGap,stdTolerance)
%% Settings
% timeGap=100; %Time in milliseconds below which individal measurements are assigned to the same slug and above which they are considered to belong to a different slug.
% stdTolerance=50; % If the relative standard deviation is above this percentage, measurements are ignored.


volumesFile=dir("*volumes.csv");
if isempty(volumesFile)
    disp("No time series volumes file found. Run imagesToVolumes first to create a time series volumes file from a folder of slug images.")
else
    T=readtable(volumesFile.name);
    disp("Analyzing file: "+volumesFile.name)

    figure
    plot(T.time, T.droplet_volume,'xk')

    [groupedTimes,groupedVolumes,T2]=groupVolumesByTime(T.time,T.droplet_volume,timeGap);

    %Only consider measurements that have a standard deviation
    %below the threshold
    indices=(T2.std_droplet_volume./T2.mean_droplet_volume*100<stdTolerance);

    T3=table(T2.time(indices),T2.mean_droplet_volume(indices),T2.std_droplet_volume(indices),'VariableNames',["time","mean_droplet_volume","std_droplet_volume"]);

    volumetricFlowRate=T3.mean_droplet_volume(2:end)/1000./minutes(diff(T3.time));
    corrTimes=T3.time(2:end);

    [V,I]=rmoutliers(volumetricFlowRate);

    T4=table(corrTimes(~I),V,'VariableNames',["time","volumetric_flow_rate"]);
    writetable(T3,volumesFile.name(1:end-4)+"_grouped.csv")
    writetable(T4,volumesFile.name(1:end-5)+"tric_flow.csv")

end
end

function [groupedTimes,groupedVolumes,results]=groupVolumesByTime(Times,Volumes,timeGap)
    indices=find((diff(Times)>milliseconds(timeGap))==1);
    [meanVolumes, stdVolumes] = deal(zeros(length(indices)-1,1));
    firstTime=duration(strings(length(indices)-1,1),'Format','hh:mm:ss.SSS');
    for j=1:length(indices)-1
        groupedTimes{j}=Times(indices(j):indices(j+1)-1);     
        groupedVolumes{j}=Volumes(indices(j):indices(j+1)-1);     
        meanVolumes(j)=mean(groupedVolumes{j});
        stdVolumes(j)=std(groupedVolumes{j});
        firstTime(j)=Times(indices(j));
        results=table(firstTime,meanVolumes,stdVolumes,'VariableNames',["time","mean_droplet_volume","std_droplet_volume"]);
    end
end