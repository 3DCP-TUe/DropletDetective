clear all; close all; clc

%% File locations and settings

%Location and filename of the load cell data
% filePath = "D:\OneDrive - TU Eindhoven\VENI - Digital Fabrication with Concrete\04_Experiments\20240927_Chapter9_1\dropletdetective\Load Cell";
% fileName = "20240930_Chapter9_1"; 
% filePath = "D:\OneDrive - TU Eindhoven\VENI - Digital Fabrication with Concrete\04_Experiments\20241007_Chapter9_2\dropletdetective\Load Cell";
% fileName = "20241007_Chapter9_2"; 
filePath = "D:\OneDrive - TU Eindhoven\VENI - Digital Fabrication with Concrete\04_Experiments\20241014_Chapter9_3\dropletdetective\Load Cell";
fileName = "20241014_Chapter9_3"; 

fileExtension = ".csv";

% Select the logger type 
% "py" refers to the OPC-UA logger in python https://github.com/arjendeetman/Python-OPC-UA 
% "ua" refers to the logger in UA Expert
loggerType = "py"; 

% Set the nozzle diameter, used to calculate yield stress from droplet mass
nozzleDiameter = 25; %mm 

% Set the bucket mass. This is a lower bound below values are not considered 
% for analysis
bucketWeight = 12; %N 

% Choose if the results should be plotted
plotting = true;

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
        meanTimeBucket(k1)=mean([firstTime2{k2}(end); firstTime2{k2}(1)]);
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
plot(meanTimeBucket,massFlowBucket)
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


%% Plot results
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
    plot(firstDur2a,diff(meanLoad),'xk')
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

saveData2=table(meanTimeBucket',massFlowBucket');
writetable(saveData2,fileName+"_ResultsMassFlowPerBucket.csv")
