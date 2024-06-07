clear all; close all; clc

%% Location and settings

%Location and filename of the load cell data
filePath = "D:\OneDrive - TU Eindhoven\VENI - Digital Fabrication with Concrete\04_Experiments\20240522_Preparation RILEM round4\Load cell\1";
fileName = "DataLogger_1"; 
fileExtension = ".csv";

%Set values for stable plateau detection
%A moving standard deviation is applied with window size "k1". When the
%value of this moving standard deviation is larger than "lim1" that
%measured value is labelled as unstable.

k1 = 10; %Window size
lim1 = 0.05; %N

interTime = 0.2; %seconds - When the time inbetween two stable measurments is larger than "interTime", the stable measurements belong to a different stable plateau.
minNo = 5; %Minimum number of stable measurements to select the stable plateau for further processing



%% Read file

cd(filePath)
T=readtable(fileName+fileExtension,'Delimiter',',');
for j=1:height(T)
    T2(j,1)=datetime(T.SourceTimeStamp{j},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','UTC');
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
    if Dur2(i)-Dur2(i-1)>seconds(interTime)
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

%TODO improve this methods. Instead of simply removing outliers, it could 
% be based on thresholds. In practice it does, however, seem to work well 
% since the outliers are located around bucket changes.
[SlugMass, I]=rmoutliers(diff(meanSlug)'); 
firstDur2a=firstDur(1:end-1);
firstTime2a=firstTime(1:end-1);
ContactTimeCorr=firstDur2a(~I)';
ContactTime=firstTime2a(~I)';


%%
fig=figure;
fig.Units='pixels';
fig.Position=[1 500 1920 500];
raw=plot(Dur,T.Value);
hold on
% xlim([minutes(1),minutes(2)])
stable=plot(Dur2,Val2,'.r');
firstDurPlot=plot(firstDur,meanSlug,'xg');
lastDurPlot=plot(lastDur,meanSlug,'xb');

fig=figure;
fig.Units='pixels';
fig.Position=[1 250 1920 500];
plot(firstDur2a,diff(meanSlug),'xk')
hold on
plot(ContactTimeCorr,SlugMass,'xr')
ylabel('Mass [N]')
disp("Outliers removed: "+sum(I)+" of "+length(SlugMass))

fig=figure;
fig.Units='pixels';
fig.Position=[1 1 1920 500];
hold on
plot(ContactTimeCorr,SlugMass,'xr')
ylabel('Mass [N]')



%% Save as csv file
saveData=table(ContactTime,ContactTimeCorr,SlugMass);
writetable(saveData,fileName+"_Results.csv")

