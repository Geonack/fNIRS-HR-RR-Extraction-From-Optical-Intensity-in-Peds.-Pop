clc;
clear;
close all;

% Choose Directory for Subject Data
cd %CD subject folder directory
allsubject= dir('*.snirf');
subjects = 1:107; 
peakcutoff = .99;
count =1;

% Peak Extraction for Each Window of Each Subject
for windows = 100:100:24000
    disp("Window: " +windows)
    alldata = [];

    % Load Raw Data
    parfor subjnumber=subjects     
        disp("Analyzing subject number: "+subjnumber)
        if subjnumber<10
            allsubjectidx=subjnumber;
        elseif subjnumber>99
            allsubjectidx=subjnumber-90;
        else
            allsubjectidx=subjnumber+8;
        end

        subjname=allsubject(allsubjectidx).name;
        raw=SnirfLoad(subjname);

        if subjnumber<17
            agedata(subjnumber, 1)=6;
        elseif subjnumber<29 && subjnumber>16
            agedata(subjnumber, 1)=7;
        elseif subjnumber<47 && subjnumber>28
            agedata(subjnumber, 1)=8;
        elseif subjnumber<59 && subjnumber>46
            agedata(subjnumber, 1)=9;
        elseif subjnumber<72 && subjnumber>58
            agedata(subjnumber, 1)=10;
        elseif subjnumber<88 && subjnumber>71
            agedata(subjnumber, 1)=11;
        elseif subjnumber<102 && subjnumber>87
            agedata(subjnumber, 1)=12;
        else
            agedata(subjnumber, 1)=13;
        end


        % Filterdata
        [data_ycRR, ~, filterout]=hmrR_PCAFilter(raw.data, [], [], 2);
        [data_ycHR, ~, filterout]=hmrR_PCAFilter(raw.data, [], [], 1);
        data_ycHR = hmrR_BandpassFilt(data_ycHR, 1, 2);
        data_ycRR = hmrR_BandpassFilt(data_ycRR, .2, .5);
        dtrendHR=detrend(data_ycHR.dataTimeSeries(3000:27000,:)); 
        dtrendRR=detrend(data_ycRR.dataTimeSeries(3000:27000,:)); 

        allPeakLocHR = [];
        allPeakLocRR = [];
        for ii = 1:92
            [pwelHR,fwelHR] = pwelch(dtrendHR(:,ii), windows, [], [], 50, 'onesided');
            [pwelRR,fRelRR] = pwelch(dtrendRR(:,ii), windows, [], [], 50, 'onesided');
            [HRpks,HRlocs] = findpeaks(pwelHR, 'MinPeakHeight', peakcutoff*max(pwelHR));
            [RRpks,RRlocs] = findpeaks(pwelRR, 'MinPeakHeight', peakcutoff*max(pwelRR));
            if length(HRlocs)~=1
                allPeakLocHR = [allPeakLocHR 0];
            else
                allPeakLocHR = [allPeakLocHR fwelHR(HRlocs)'];
            end
            if length(RRlocs)~=1
                allPeakLocRR = [allPeakLocRR 0];
            else
                allPeakLocRR = [allPeakLocRR fRelRR(RRlocs)'];
            end
        end
        allDataHR(:,subjnumber,count)=allPeakLocHR';
        allDataRR(:,subjnumber,count)=allPeakLocRR';
    end

    for kk = 1:107
        HRanalyzisdata = allDataHR(allDataHR(:,kk)~=0,kk);
        RRanalyzisdata = allDataRR(allDataRR(:,kk)~=0,kk);
        HRout = isoutlier(HRanalyzisdata,"mean");
        RRout = isoutlier(RRanalyzisdata,"mean");
        HRanalyzisdata = HRanalyzisdata(~HRout);
        RRanalyzisdata = RRanalyzisdata(~RRout);
        HRmeans(kk) = mean(HRanalyzisdata);
        RRmeans(kk) = mean(RRanalyzisdata);
    end
    % Count Data For Graphs
    n(count,1) = windows;
    n(count,2) = length(unique(HRmeans));
    n(count,3) = length(unique(RRmeans));
    count = count +1;
end

save("ChanByChanHRCorr.mat")

%% Correlation values for each window HR
for window = 1:240
    disp("Window: " + window.*100)
    for kk = 1:107
        analyzisdata = allDataHR(allDataHR(:,kk,window)~=0,kk,window);
        out = isoutlier(analyzisdata,"mean");
        analyzisdata = analyzisdata(~out);
        means(kk) = mean(analyzisdata);
    end

    for ll = 1:107
        disp("Analyzing subject number: "+ll)
        if ll<10
            allsubjectidx=ll;
        elseif ll>99
            allsubjectidx=ll-90;
        else
            allsubjectidx=ll+8;
        end

        subjname=allsubject(allsubjectidx).name;
        raw=SnirfLoad(subjname);
        rawdata = hmrR_BandpassFilt(raw.data,1,2);
        meanband = hmrR_BandpassFilt(raw.data,means(ll),means(ll));
        for mm = 1:92
            meancorr(mm,ll) = corr(rawdata.dataTimeSeries(3000:27000,mm),meanband.dataTimeSeries(3000:27000,mm));
        end
    end
    allmeansHR(window,1) = window*100;
    allmeansHR(window,2) = mean(meancorr,"all");
    allmeansHR(window,3) = median(meancorr,"all");
    allmeansHR(window,4) = mode(meancorr,"all");
    allmeansHR(window,5) = range(meancorr,"all");
end



%% Correlation values each window RR

for window = 1:240
    disp("Window: " + window.*100)
    for kk = 1:107
        analyzisdata = allDataRR(allDataRR(:,kk,window)~=0,kk,window);
        out = isoutlier(analyzisdata,"mean");
        analyzisdata = analyzisdata(~out);
        means(kk) = mean(analyzisdata);
    end

    for ll = 1:107
        disp("Analyzing subject number: "+ll)
        if ll<10
            allsubjectidx=ll;
        elseif ll>99
            allsubjectidx=ll-90;
        else
            allsubjectidx=ll+8;
        end

        subjname=allsubject(allsubjectidx).name;
        raw=SnirfLoad(subjname);
        rawdata = hmrR_BandpassFilt(raw.data,.2,.5);
        meanband = hmrR_BandpassFilt(raw.data,means(ll),means(ll));
        for mm = 1:92
            meancorr(mm,ll,window) = corr(rawdata.dataTimeSeries(3000:27000,mm),meanband.dataTimeSeries(3000:27000,mm));
        end
    end
    allmeansRR(window,1) = window*100;
    allmeansRR(window,2) = mean(meancorr,"all");
    allmeansRR(window,3) = median(meancorr,"all");
    allmeansRR(window,4) = mode(meancorr,"all");
    allmeansRR(window,5) = range(meancorr,"all");
end

%% Graphing

regres = "fourier8";
for window = 1:240
    HRests = allDataHR(:,:,window);
    HRincl = HRests(HRests~=0);
    allwindow(window) = window*100;
    excluHR(window) = ((107*92)-length(HRincl))./(107*92);
    uniqvalsHR(window)= length(unique(allDataHR(:,:,window)))
end
fi = fit(allwindow',excluHR',regres);
fitlength = linspace(2^14,24000,10000)
fitdata = feval(fi,fitlength);
[values,troughs] = findpeaks(-fitdata,'SortStr','descend','NPeaks',1);
for ii = 1:round(feval(fi,fitlength(troughs))*10000)
    line(ii,1) = round(fitlength(troughs),-2);
    line(ii,2) = ii/10000;
end

%***********************************************************************************************************************
hfig = figure;

% Create axes
axes1 = axes('Parent',hfig);
hold(axes1,'on');

% Create plot
plot(allwindow,excluHR,'DisplayName','data','MarkerSize',20,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 1]);

% Create plot
plot(allwindow,feval(fi,allwindow),'DisplayName','fitted curve','LineWidth',5,'Color',[1 0 0]);

% Create plot
plot(line(:,1),line(:,2),'LineWidth',5,'LineStyle','--','Color',[0 0 0]);

% Create ylabel
ylabel('Percentage of Channels Excluded Across All Subjects','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create xlabel
xlabel('Window Size (Samples)','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create title
title('Heart Rate','FontWeight','bold','FontName','Times New Roman','FontSize',20);

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 24000]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 .06]);
box(axes1,'on');
hold(axes1,'off');


xticklabels(["0","5,000","10,000","15,000","20,000"])
axis("square")
%%
hfig = figure;

% Create axes
axes1 = axes('Parent',hfig);
hold(axes1,'on');

% Create plot
plot(allwindow',uniqvalsHR,'LineWidth',5,'Color',[0 0 0]);

% Create ylabel
ylabel('Number of unique Values','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create xlabel
xlabel('Window Size (Samples)','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create title
title('Heart Rate','FontName','Times New Roman','FontSize',20);

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
ylim([0 600])
xlim([0 24000])

xticklabels(["0","5,000","10,000","15,000","20,000"])
axis("square")

%%
for window = 1:240
    RRests = allDataRR(:,:,window);
    RRincl = RRests(RRests~=0);
    allwindow(window) = window*100
    excluRR(window) = ((107*92)-length(RRincl))./(107*92)
    uniqvalsRR(window)= length(unique(allDataRR(:,:,window)))
end
fi2 = fit(allwindow',excluRR',regres);
fitlength2 = linspace(2^14,24000,10000)
fitdata2 = feval(fi2,fitlength2);
[values2,troughs2] = findpeaks(-fitdata2,'SortStr','descend','NPeaks',1);
for ii = 1:round(feval(fi2,fitlength2(troughs2)).*10000)
    line2(ii,1) = round(fitlength2(troughs2),-2);
    line2(ii,2) = ii/10000;
end

%***********************************************************************************************************************
hfig = figure;

% Create axes
axes1 = axes('Parent',hfig);
hold(axes1,'on');

% Create plot
plot(allwindow,excluRR,'DisplayName','data','MarkerSize',20,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 1]);

% Create plot
plot(allwindow,feval(fi2,allwindow),'DisplayName','fitted curve','LineWidth',5,'Color',[1 0 0]);

% Create plot
plot(line2(:,1),line2(:,2),'LineWidth',5,'LineStyle','--','Color',[0 0 0]);

% Create ylabel
ylabel('Percentage of Channels Excluded Across All Subjects','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create xlabel
xlabel('Window Size (Samples)','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create title
title('Respiration Rate','FontWeight','bold','FontName','Times New Roman','FontSize',20);

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 24000]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 .06]);
box(axes1,'on');
hold(axes1,'off');

xticklabels(["0","5,000","10,000","15,000","20,000"])
axis("square")

%%
hfig = figure;

% Create axes
axes1 = axes('Parent',hfig);
hold(axes1,'on');

% Create plot
plot(allwindow',uniqvalsRR,'LineWidth',5,'Color',[0 0 0]);

% Create ylabel
ylabel('Number of unique Values','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create xlabel
xlabel('Window Size (Samples)','FontWeight','bold',...
    'FontName','Times New Roman','FontSize',20);

% Create title
title('Respiration Rate','FontName','Times New Roman','FontSize',20);

box(axes1,'on');
hold(axes1,'off');
xlim(axes1,[0 24000]);
ylim([0 600])
xticklabels(["0","5,000","10,000","15,000","20,000"])
axis("square")