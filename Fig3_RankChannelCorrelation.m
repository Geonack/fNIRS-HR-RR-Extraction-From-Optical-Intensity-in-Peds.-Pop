clc;
clear;
close all;

% Choose Directory for Subject Data
cd %CD subject folder directory
allsubject= dir('*.snirf');
subjects = 1:107;
%% Load Raw Data
for subjnumber=subjects     
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

    %bandpass filter data
    stdata = hmrR_BandpassFilt(raw.data,.2,2);
    for i = 1:92
        stdev(i,1) = std(diff(zscore(detrend(stdata.dataTimeSeries(:,i)))));
        stdev(i,2) = i;
    end
    Spearman(:,subjnumber)=stdev(:,1);
    sort1 = sortrows(stdev,1);
    AllChanRank(:,subjnumber) = sort1(:,2);

    %Rank Channels
    for i = 1:46
        stdev2(i,1) = std(diff(zscore(detrend(stdata.dataTimeSeries(:,i)))));
        stdev2(i+46,1) = std(diff(zscore(detrend(stdata.dataTimeSeries(:,i+46)))));
        stdev2(i,2) = i;
        stdev2(i+46,2) = i;
    end
    sort2_690 = sortrows(stdev2(1:46,:),1);
    sort2_830 = sortrows(stdev2(47:end,:),1);
    for ii = 1:length(sort2_690)
        listsort_690(ii) = find(sort2_690 == ii)-46;
        listsort_830(ii) = find(sort2_830 == ii)-46;
    end

    AllChanRank690(:,subjnumber) = listsort_690';
    AllChanRank830(:,subjnumber) = listsort_830';

end

%% Final figure 1

%Define left and right hemisphere channels
leftHem= [1 2 3 4 9 10 11 12 13 14 15 16 25 26 27 28 29 30 31 39 40 41 42];
rightHem= [5 6 7 8 17 18 19 20 21 22 23 24 32 33 34 35 36 37 38 43 44 45 46];

%Left Hem 690
hfig = figure;
t = 1:46;
textsize = 20;
avg = mean(AllChanRank690');
med = median(AllChanRank690');
polylinewidth = 10;
avglinewidth = 5;

plot(t(leftHem),avg(leftHem),"o","MarkerFaceColor","r","MarkerEdgeColor","k","MarkerSize",polylinewidth)
hold on
plot(t(leftHem),med(leftHem),"^","MarkerFaceColor","b","MarkerEdgeColor","k","MarkerSize",polylinewidth)
ff = fit(t(leftHem)',avg(leftHem)',"poly1");
plots = plot(ff,"r:");
plots.LineWidth = avglinewidth;
ff = fit(t(leftHem)',med(leftHem)',"poly1");
plots = plot(ff,"b--");
plots.LineWidth = avglinewidth;
hold off
xlabel("Channel")
ylabel("Rank")
xlim([0 48]);
ylim([0 48]);
legend("Mean", "Median","Mean","Median","Location","eastoutside")

[R690LA,P690LA] = corrcoef(t(leftHem),avg(leftHem))
[R690LM,P690LM] = corrcoef(t(leftHem),med(leftHem))

title("690 nm Left Hemisphere")
axis square
set(findall(hfig,'-property','FontSize'),'FontSize',textsize) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%%Right Hem 690
hfig = figure;
t = 1:46;
textsize = 20;
avg = mean(AllChanRank690');
med = median(AllChanRank690');
polylinewidth = 10;
avglinewidth = 5;

plot(t(rightHem),avg(rightHem),"o","MarkerFaceColor","r","MarkerEdgeColor","k","MarkerSize",polylinewidth)
hold on
plot(t(rightHem),med(rightHem),"^","MarkerFaceColor","b","MarkerEdgeColor","k","MarkerSize",polylinewidth)
ff = fit(t(rightHem)',avg(rightHem)',"poly1");
plots = plot(ff,"r:");
plots.LineWidth = avglinewidth;
ff = fit(t(rightHem)',med(rightHem)',"poly1");
plots = plot(ff,"b--");
plots.LineWidth = avglinewidth;
hold off
xlabel("Channel")
ylabel("Rank")
xlim([0 48]);
ylim([0 48]);
legend("Mean", "Median","Mean","Median","Location","eastoutside")

[R690RA,P690RA] = corrcoef(t(rightHem),avg(rightHem))
[R690RM,P690RM] = corrcoef(t(rightHem),med(rightHem))

title("690 nm Right Hemisphere")

axis square
set(findall(hfig,'-property','FontSize'),'FontSize',textsize) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%% Final figure 2
%Heatmap of all channels based on approximate position
channelMap = [0 0 1 2 0 5 6 0 0 0 0 3 4 0 7 8 0 0 9 10 13 14 0 17 18 21 22 11 12 15 16 0 19 20 23 24 ...
    25 26 28 29 0 32 33 36 37 0 27 30 31 0 34 35 38 0 0 0 29 40 0 43 44 0 0 0 0 41 42 0 45 46 0 0 ];
dispMap = NaN([1 72]);

for i = 1:length(channelMap)
    if channelMap(i) == 0
    else 
        dispMap(i) =avg(channelMap(i))
    end
end

dispMap = reshape(dispMap,[9 8]);
dispMap = dispMap';

figure;
h = heatmap(dispMap,"FontSize",5,"FontName","Times")
h.Colormap = turbo(256)
clim([14 31])
h.MissingDataColor = [1 1 1]
grid off

% 830 NM
%Left Hem 830
hfig = figure;
t = 1:46;
textsize = 20;
avg = mean(AllChanRank830');
med = median(AllChanRank830');
polylinewidth = 10;
avglinewidth = 5;

plot(t(leftHem),avg(leftHem),"o","MarkerFaceColor","r","MarkerEdgeColor","k","MarkerSize",polylinewidth)
hold on
plot(t(leftHem),med(leftHem),"^","MarkerFaceColor","b","MarkerEdgeColor","k","MarkerSize",polylinewidth)
ff = fit(t(leftHem)',avg(leftHem)',"poly1");
plots = plot(ff,"r:");
plots.LineWidth = avglinewidth;
ff = fit(t(leftHem)',med(leftHem)',"poly1");
plots = plot(ff,"b--");
plots.LineWidth = avglinewidth;
hold off
xlabel("Channel")
ylabel("Rank")
xlim([0 48]);
ylim([0 48]);
legend("Mean", "Median","Mean","Median","Location","eastoutside")

[R830LA,P690LA] = corrcoef(t(leftHem),avg(leftHem))
[R830LM,P830LM] = corrcoef(t(leftHem),med(leftHem))

title("830 nm Left Hemisphere")

axis square
set(findall(hfig,'-property','FontSize'),'FontSize',textsize) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%%Right Hem 830
hfig = figure;
t = 1:46;
textsize = 20;
avg = mean(AllChanRank830');
med = median(AllChanRank830');
polylinewidth = 10;
avglinewidth = 5;

plot(t(rightHem),avg(rightHem),"o","MarkerFaceColor","r","MarkerEdgeColor","k","MarkerSize",polylinewidth)
hold on
plot(t(rightHem),med(rightHem),"^","MarkerFaceColor","b","MarkerEdgeColor","k","MarkerSize",polylinewidth)
ff = fit(t(rightHem)',avg(rightHem)',"poly1");
plots = plot(ff,"r:");
plots.LineWidth = avglinewidth;
ff = fit(t(rightHem)',med(rightHem)',"poly1");
plots = plot(ff,"b--");
plots.LineWidth = avglinewidth;
hold off
xlabel("Channel")
ylabel("Rank")
xlim([0 48]);
ylim([0 48]);
legend("Mean", "Median","Mean","Median","Location","eastoutside")

[R830RA,P830RA] = corrcoef(t(rightHem),avg(rightHem))
[R830RM,P830RM] = corrcoef(t(rightHem),med(rightHem))

title("830 nm Right Hemisphere")
%ylim([0 46]);

axis square
set(findall(hfig,'-property','FontSize'),'FontSize',textsize) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%830 mn
for i = 1:length(channelMap)
    if channelMap(i) == 0
    else 
        dispMap(i) =avg(channelMap(i))
    end
end

dispMap = reshape(dispMap,[9 8]);
dispMap = dispMap';

figure;
h = heatmap(dispMap,"FontSize",20,"FontName","Times")
h.Colormap = turbo(256)
clim([14 31])
h.MissingDataColor = [1 1 1]
grid off

%% Final figure 3 Heatmap of Channel comparisons based on correlation rank
textsize = 20;
hfig = figure;
r = corr(Spearman(1:46,:)',"type","Spearman");
clims = [-1 1];
imagesc(r,clims)
a = colorbar;
colormap(redblue(256));
ylabel(a,"Spearman's rho",'Interpreter','latex');
xlabel("Channels")
ylabel("Channels")
title("690 nm")


set(findall(hfig,'-property','FontSize'),'FontSize',textsize)
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
axis("square")
hfig = figure;
r2 = corr(Spearman(47:end,:)',"type","Spearman");
clims = [-1 1];
imagesc(r2,clims)
a = colorbar;
colormap(redblue(256));
ylabel(a,"Spearman's rho",'Interpreter','latex');
xlabel("Channels")
ylabel("Channels")
title("830 nm")
set(findall(hfig,'-property','FontSize'),'FontSize',textsize)
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
axis("square")

