clear all; close all; clc

%% Location and settings

%Location and filename of the load cell data
filePath = "";
fileName = ""; 
fileExtension = ".csv";
loggerType = "py"; %"py" refers to the OPC-UA logger in python https://github.com/arjendeetman/Python-OPC-UA or "ua" refers to the logger in UA Expert

nozzleDiameter = 25; %mm

%Set values for stable plateau detection
%A moving standard deviation is applied with window size "k1". When the
%value of this moving standard deviation is larger than "lim1" that
%measured value is labelled as unstable.

k1 = 10; %Window size
lim1 = 0.05; %N

interTimeLimit = 0.2; %seconds - When the time inbetween two stable measurments is larger than "interTime", the stable measurements belong to a different stable plateau.
minNo = 5; %Minimum number of stable measurements to select the stable plateau for further processing

filterMethod = "medianDetection"; %removeOutliers or "medianDetection"
%Setting for filterMethod median detection
intervalLimit=0.7;

plotting = true;

%The limits below are only used for filterMethod "setLimits"


%% Read file

cd(filePath)
T=readtable(fileName+fileExtension,'Delimiter',',');
if loggerType=="ua"
    for j=1:height(T)
        T2(j,1)=datetime(T.SourceTimeStamp{j},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','UTC');
    end
elseif loggerType=="py"
    T2=T.Time;
    T.Value=T.Load;
end

%%
Dur=T2-T2(1); 
S=movstd(T.Value,k1);

Val2=T.Value(S<lim1); %Stable measurements
Dur2=Dur(S<lim1); %Duration corresponding to stable measurements
Time2=T2(S<lim1); %Times corresponding to stable measurements

clear S Val3 Dur3 firstDur lastDur meanSlug firstTime lastTime
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
            meanSlug(k)=mean(Val3{k});
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

interTime=firstTime(2:end)-lastTime(1:end-1);

% TODO improve this methods. Instead of simply removing outliers, it could 
% be based on thresholds. In practice it does, however, seem to work well 
% since the outliers are located around bucket changes.
if filterMethod =="removeOutliers"
    [SlugMass, I]=rmoutliers(diff(meanSlug)');
    firstDur2a=firstDur(2:end);
    firstTime2a=firstTime(2:end);
    ContactTimeCorr=firstDur2a(~I)';
    ContactTime=firstTime2a(~I)';
    disp("Outliers removed: "+sum(I)+" of "+length(SlugMass))
elseif filterMethod == "medianDetection"
    SlugMassA=diff(meanSlug)';
    SlugMassB=SlugMassA(interTime>(1-intervalLimit)*median(interTime)&interTime<(1+intervalLimit)*median(interTime));
    SlugMass=SlugMassB(SlugMassB>0.1*median(SlugMassB));
    % SlugMass=SlugMassA(interTime>0.3*mode(SlugMassA)&SlugMassA<1.7*mode(SlugMassA));

    firstDur2a=firstDur(2:end);
    firstTime2a=firstTime(2:end);
    ContactTimeCorrA=firstDur2a(interTime>(1-intervalLimit)*median(interTime)&interTime<(1+intervalLimit)*median(interTime));
    ContactTimeCorr=ContactTimeCorrA(SlugMassB>0.1*median(SlugMassB))';

    ContactTimeA=firstTime2a(interTime>(1-intervalLimit)*median(interTime)&interTime<(1+intervalLimit)*median(interTime));
    ContactTime=ContactTimeA(SlugMassB>0.1*median(SlugMassB))';
end
YieldStress=SlugMass/(sqrt(3)*1/4*pi*nozzleDiameter^2)*1000; %kPa
disp("Sum "+sum(SlugMass)/9.81+" kg - "+length(SlugMass) +" slugs")
 
% figure
% plot(SlugMassA,'xk')
% hold on
% plot(SlugMassB,'xr')


%%
if plotting == true
    close all

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 500 1920 500];
    raw=plot(Dur,T.Value);
    hold on
    % xlim([minutes(1),minutes(2)])
    stable=plot(Dur2,Val2,'.r');
    % firstDurPlot=plot(firstDur,meanSlug,'xg');
    % lastDurPlot=plot(lastDur,meanSlug,'xb');
    yyaxis right
    hold on
    pl=plot(ContactTimeCorr,SlugMass,'x','Color','m');
    ax=gca;
    ax.YAxis(2).Color='m';

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 250 1920 500];
    grid on
    box on
    plot(firstDur2a,diff(meanSlug),'xk')
    hold on
    plot(ContactTimeCorr,SlugMass,'xr')
    ylabel('Mass [N]')

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 1 1920 500];
    hold on
    plot(ContactTimeCorr,YieldStress,'xr')
    ylabel('Yield stress [kPa]')

    fig=figure;
    fig.Units='pixels';
    fig.Position=[1 1 1920/3 500];
    hold on
    histogram(YieldStress,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5)
    xlabel('Yield stress [kPa]')
    ylabel('Frequency')
end

%% Save as csv file
saveData=table(ContactTime,ContactTimeCorr,SlugMass,YieldStress);
writetable(saveData,fileName+"_Results.csv")

