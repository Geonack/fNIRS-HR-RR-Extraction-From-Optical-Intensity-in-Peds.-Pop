clc;
clear;
close all;


%Load Subject folder
cd %CD subject folder directory
allsubject= dir('*.snirf');
subjnumber = 107; %Change to change subject extracted
sampchannel = 23; %Change to change sample extracted

if sampchannel >=46
    chandisp = sampchannel-46;
    wavelength = 830;
else
    chandisp = sampchannel;
    wavelength = 690;
end

if subjnumber<10
    allsubjectidx=subjnumber;
elseif subjnumber>99
    allsubjectidx=subjnumber-90;
else
    allsubjectidx=subjnumber+8;
end

%Extract Raw Data
subjname=allsubject(allsubjectidx).name;
raw=SnirfLoad(subjname);
   
intervalStart = 11000;%Select interval of time to start displaying from
intervalLength = 30; %Select interval of time in seconds to display
interval = intervalStart:intervalStart+(intervalLength/0.02);

rawRaw = raw.data;
times = 0:0.02:intervalLength;
[data_pcHR, ~, filterout]=hmrR_PCAFilter(raw.data, [], [], 1);
[data_pcRR, ~, filterout]=hmrR_PCAFilter(raw.data, [], [], 2);
data_ycHR = hmrR_BandpassFilt(data_pcHR, 1, 2);
data_ycRR = hmrR_BandpassFilt(data_pcRR, .2, .5);
dtrendHR=detrend(data_ycHR.dataTimeSeries(3000:27000,:)); 
dtrendRR=detrend(data_ycRR.dataTimeSeries(3000:27000,:));
[pwelHR,fwelHR] = pwelch(dtrendHR(:,sampchannel), 17700, [], [], 50, 'onesided');
[pwelRR,fwelRR] = pwelch(dtrendRR(:,sampchannel), 17400, [], [], 50, 'onesided');
[HRpks,HRlocs] = findpeaks(pwelHR, 'MinPeakHeight', .99*max(pwelHR));
[RRpks,RRlocs] = findpeaks(pwelRR, 'MinPeakHeight', .99*max(pwelRR));
fwelHR(HRlocs)
fwelRR(RRlocs)

%% Raw Data Figure
hfig = figure;
sizefont = 15;
sizefont2 = 12;
plot(times,rawRaw.dataTimeSeries(interval,sampchannel),"k")
xlim([times(1) times(end)])
xlabel("Time (sec)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel('Optical Intensity (W/m^2)',"FontSize",sizefont2,"FontName","Times","FontWeight","bold")
picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'rawdatafig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)
%% PCA HR Figure
hfig = figure;
plot(times,data_pcHR.dataTimeSeries(interval,sampchannel),"k")
xlim([times(1) times(end)])
xlabel("Time (sec)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("Optical Intensity (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'PCAHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)

%% PCA RR Figure 
hfig = figure;
plot(times,data_pcRR.dataTimeSeries(interval,sampchannel),"k")
xlim([times(1) times(end)])
xlabel("Time (sec)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("Optical Intensity (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")

picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'PCAHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)
%% DETREND HR Figure
hfig = figure;
plot(times,dtrendHR(interval-3000,sampchannel),"k","LineWidth",1.25)
xlim([times(1) times(end)])
xlabel("Time (sec)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("Optical Intensity (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")

picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'detrendHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)

%% DETREND RR Figure
hfig = figure;
plot(times,dtrendRR(interval-3000,sampchannel),"k","LineWidth",1.25)
xlim([times(1) times(end)])
xlabel("Time (sec)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("Optical Intensity (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")

picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'detrendHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)


%% PWelch HR Figure
hfig = figure;
plot(fwelHR,pwelHR,"k","LineWidth",1.25)
hold on
plot(fwelHR(HRlocs),HRpks,"r.","MarkerSize",25)
hold off
xlim([0 3])
xlabel("Frequency (Hz)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("PSD (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")

picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'psdHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)

%% Plot RR
hfig = figure;
plot(fwelRR,pwelRR,"k","LineWidth",1.25)
hold on
plot(fwelRR(RRlocs),RRpks,"r.","MarkerSize",25)
hold off
xlim([0 3])
xlabel("Frequency (Hz)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")
ylabel("PSD (A.U.)","FontSize",sizefont2,"FontName","Times","FontWeight","bold")

picturewidth = 15;
set(findall(hfig,'-property','Box'),'Box','off')
set(hfig,'Units','centimeters','Position',[1 1 picturewidth 7.5])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig,'psdHRfig.pdf','Resolution',300,'ContentType','vector')
set(gca,'fontsize',14)