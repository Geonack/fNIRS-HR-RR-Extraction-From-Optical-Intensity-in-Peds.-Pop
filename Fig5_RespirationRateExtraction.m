clc;
clear;
close all;

%Load Demographic Data Folder
cd %CD demographic folder directory
x = readtable("ChildDatatset.csv");
sex = x{2:end,6};
sex = int2str(sex)
sex(sex == '1') ="M"
sex(sex == '2') ="F"
nominalAge = x{2:end,9};

%Load fNIRS data folder
cd %CD subject folder directory
allsubject= dir('*.snirf');
subjects = 1:107;
peakcutoff = .99;
prca = 2;
%% Load Raw Data
for subjnumber= subjects     
    % Create Subject Number Index
    disp("Analyzing subject number: "+subjnumber)
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

    %Obtain Age Data
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
    
    %% Filterdata

   
    [data_ycHR, ~, filterout]=hmrR_PCAFilter(raw.data, [], [], prca);
    data_ycHR = hmrR_BandpassFilt(data_ycHR, .2, .5);
    dtrendHR=detrend(data_ycHR.dataTimeSeries(3000:27000,:)); 
    
    %Perform Welch's method and optain peak locations
    allPeakLoc = [];
    parfor ii = 1:92
        [pwel,fwel] = pwelch(dtrendHR(:,ii),6400 , [], [], 50, 'onesided');
        [welpks,wellocs] = findpeaks(pwel, 'MinPeakHeight', peakcutoff*max(pwel));
        if length(wellocs)~=1
            allPeakLoc = [allPeakLoc 0];
        else
            allPeakLoc = [allPeakLoc fwel(wellocs)'];
        end
    end
    allData(:,subjnumber)=allPeakLoc';
end

%% Remove Outliers
parfor kk = 1:107
    analyzisdata = allData(allData(:,kk)~=0,kk);
    out = isoutlier(analyzisdata,"mean");
    analyzisdata = analyzisdata(~out);
    totalChansExcluded(kk) = 92-length(analyzisdata);
    means(kk) = mean(analyzisdata);
    stanDevs(kk) = std(analyzisdata)
end

%% HR vs Age Figure

hfig = figure();
p = plot(nominalAge(sex == 'M'),means(sex == 'M'),"o","MarkerFaceColor","r","MarkerEdgeColor","k");
p.MarkerSize = 7;
hold on
p = plot(nominalAge(sex == 'F'),means(sex == 'F'),"^","MarkerFaceColor","b","MarkerEdgeColor","k");
p.MarkerSize = 7;

ff = fit(nominalAge(sex == 'M'),means(sex == 'M')',"poly1");
[R1,P1] = corrcoef(nominalAge(sex == 'M'),means(sex == 'M')');
disp("Mean R: "+string(R1(1,2))+" P: " +string(P1(1,2)))
hold on
p = plot(ff,"r--");
p.LineWidth = 4;

ff = fit(nominalAge(sex == 'F'),means(sex == 'F')',"poly1");
[R1,P1] = corrcoef(nominalAge(sex == 'F'),means(sex == 'F')');
disp("Mean R: "+string(R1(1,2))+" P: " +string(P1(1,2)))
hold on
p = plot(ff,"b:");
p.LineWidth = 4;
hold off

ff = fit(nominalAge,means',"poly1");
[R1,P1] = corrcoef(nominalAge,means');
disp("Mean R: "+string(R1(1,2))+" P: " +string(P1(1,2)))
hold on
p = plot(ff,"k");
p.LineWidth = 4;


xlim([5.5 14.5]);
xlabel("Age (years)",'FontSize',20,"FontName","Times",'FontWeight','bold')
ylabel("Respiration Rate (Hz)",'FontSize',20,"FontName","Times",'FontWeight','bold')
legend("Male Data", "Female Data","Male Data Regression", ...
    "Female Data Regression","All Data Regression", ...
    'FontSize',15,"Location","eastoutside","FontName","Times",'FontWeight','bold')

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)

%% Variablity in subjects plot
parfor kk = subjects
    intData = allData(:,kk)
    intData(intData==0)=NaN;
    totalData(:,kk)=intData;
end
hfig = figure;
boundedline(1:107,means,stanDevs,'k')
hold on
plot(1:107,means,"k","LineWidth",2)
hold off
xlim([0 max(subjects)+1]);
xlabel("Subjects","FontName","Times",'FontWeight','bold')
ylabel("RR (Hz)","FontName","Times",'FontWeight','bold')

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)

%% T-test Male vs female
hfig = figure();
mData = means(sex == 'M');
fData = means(sex == 'F');

[hRes,pRes,ciRes,statRes] = ttest2(mData,fData)

histogram(mData*60,[16:1:22])
hold on
histogram(fData*60,[16:1:22])
hold off
legend("Male","Female","FontName","Times",'FontWeight','bold')

xlabel("Respiration Rate (Breaths/Min)",'FontSize',20,"FontName","Times",'FontWeight','bold')
ylabel("Frequency",'FontSize',20,"FontName","Times",'FontWeight','bold')

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)

%% Intra-group Variability Graph
hfig = figure();
boxplot(means,agedata,'Colors','k','OutlierSize',10,'Symbol','*','Widths',0.8);

xlabel("Age (years)",'FontSize',20,"FontName","Times",'FontWeight','bold')
ylabel("RR (Hz)",'FontSize',20,"FontName","Times",'FontWeight','bold')


set(findall(hfig,'-property','Box'),'Box','off') % optional
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)