clc;
clear;
close all;

%Uncomment whether HR or RR lags should be extracted
type = "HR"
%type = "RR"

%Uncomment whether 690 or 830 nm lags should be extracted
nm = 690
%nm = 830

%%Load Extracted Lags for LagExtractions.m
if type == "HR"
    load("AbsoluteValueLagsHR.mat")
else
    load("AbsoluteValueLagsRR.mat")
end

if nm == 690
    ranges = 1:46;
else
    ranges= 47:92;
end

%Get Ages
x = readtable("ChildDatatset.csv");
sex = x{2:end,6};
sex = int2str(sex)
sex(sex == '1') ="M"
sex(sex == '2') ="F"

    
plotlags = meanlags';

%% Repeated Measures Anova
close all
ages2 = 5:14;
formatSpec = 'Channel%d';
for A1 = 1:46
    channames(A1) = string(sprintf(formatSpec,A1));
end
ages = [ones(1,16).*6 ones(1,12).*7 ones(1,18).*8 ones(1,12).*9 ones(1,13).*10 ones(1,16).*11 ones(1,14).*12 ones(1,6).*13]';
lags = array2table(plotlags(:,ranges),'VariableNames',channames);
lags = addvars(lags,string(ages),'Before','Channel1','NewVariableNames',"Age");
Meas = table([1:46]',VariableNames="Channels");
rm = fitrm(lags,'Channel1-Channel46~Age','WithinDesign',Meas);
ranovatbl = ranova(rm)
spheric = mauchly(rm)
epsilons = epsilon (rm)


tbl = margmean(rm,'Age','alpha',.05);
tbl.Age=str2double(tbl{:,1});
tbl = sortrows(tbl,"Age");
xlim([5 14]);
close all
figure1 = figure();
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

bar(tbl.Age,tbl.Mean);
hold on
errorbar(tbl.Age,tbl.Mean,tbl.Lower,tbl.Upper,"LineWidth",3,"LineStyle","none");
hold off

% Create ylabel
ylabel('Delay (Seconds)','FontWeight','bold','FontName','Times New Roman','FontSize',25);

% Create xlabel
xlabel('Age (Years)','FontWeight','bold','FontName','Times New Roman','FontSize',25);


xlim(axes1,[5.5 13.5]);
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontName','Times New Roman','FontSize',15);

picturewidth = 20;
hw_ratio = 0.65;
set(findall(figure1,'-property','Box'),'Box','off') 
set(findall(figure1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figure1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(figure1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square


%% Multicompare Ages
figure;
comp = multcompare(rm,'Age','ComparisonType','lsd');
comp.Age_1=str2double(comp{:,1});
comp.Age_2=str2double(comp{:,2});
comp = sortrows(comp,"Age_1");
comp_matrix = NaN(13);
for i = 1:length(comp.Age_1)
    comp_matrix(comp{i,1},comp{i,2})=comp{i,5};
end
comp_matrix = comp_matrix(6:13,6:13);

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(comp_matrix,.05,'dep','yes');
subplot(1,2,1)
imagesc(comp_matrix);
axis square
colorbar;
title("Paired T No Correction")
subplot(1,2,2)
imagesc(h);
axis square
colorbar;
title("Significant after FDR_B_H (1 is significant)")

[row,col] = find(h==1) 
indices = sortrows([row+5 col+5],1)

%% Age comparison means chart

graphicmeans(1) = mean(lags{1:16,2:end},"all");%6
graphicmeans(2) = mean(lags{17:28,2:end},"all");%7
graphicmeans(3) = mean(lags{29:46,2:end},"all");%8
graphicmeans(4) = mean(lags{47:58,2:end},"all");%9
graphicmeans(5) = mean(lags{59:71,2:end},"all");%10
graphicmeans(6) = mean(lags{72:87,2:end},"all");%11
graphicmeans(7) = mean(lags{88:101,2:end},"all");%12
graphicmeans(8) = mean(lags{102:end,2:end},"all");%13


graphicstes(1) = std(lags{1:16,2:end},0,"all")./sqrt(46*min(size(lags{1:16,2:end})));%6
graphicstes(2) = std(lags{17:28,2:end},0,"all")./sqrt(46*min(size(lags{17:28,2:end})));%7
graphicstes(3) = std(lags{29:46,2:end},0,"all")./sqrt(46*min(size(lags{29:46,2:end})));%8
graphicstes(4) = std(lags{47:58,2:end},0,"all")./sqrt(46*min(size(lags{47:58,2:end})));%9
graphicstes(5) = std(lags{59:71,2:end},0,"all")./sqrt(46*min(size(lags{59:71,2:end})));%10
graphicstes(6) = std(lags{72:87,2:end},0,"all")./sqrt(46*min(size(lags{72:87,2:end})));%11
graphicstes(7) = std(lags{88:101,2:end},0,"all")./sqrt(46*min(size(lags{88:101,2:end})));%12
graphicstes(8) = std(lags{102:end,2:end},0,"all")./sqrt(46*min(size(lags{102:end,2:end})));%13


figure1 = figure();
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(6:13,graphicmeans,"FaceColor","white","BarWidth",1)
hold on
errorbar(6:13,graphicmeans,graphicstes,"LineStyle","None",Color="k",LineWidth=2);
if type == "HR" && nm == 690
 sig = sigstar({[6,8],[6,9],[6,10],[6,11],[6,12],[6,13],[7,9],[7,11],[7,12],[7,13]}, ...
    [0.0295 0.0040 0.0295 0.0028 0.0033 0.0093 0.0331 0.0252 0.0295 0.0391]); %HR 690 NM
end
hold off
% Create ylabel
ylabel('Delay (Seconds)','FontWeight','bold','FontName','Times New Roman','FontSize',25);
% Create xlabel
xlabel('Age (Years)','FontWeight','bold','FontName','Times New Roman','FontSize',25);
xlim(axes1,[5.5 13.5]);
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontName','Times New Roman','FontSize',15);
picturewidth = 20; 
hw_ratio = 0.65; 
set(findall(figure1,'-property','Box'),'Box','off') 
set(findall(figure1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figure1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(figure1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)
if type == "HR" && nm == 690
    title("HR 690 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 0.5])
elseif type == "HR" && nm == 830
    title("HR 830 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 0.5])
elseif type == "RR" && nm == 690
    title("RR 690 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 2.5])
else
    title("RR 830 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 2.5])
end

%% Intra-group Variability Graph
hfig = figure();

means= mean(lags{1:end,2:end},2)
boxplot(means,ages,'Colors','k','OutlierSize',10,'Symbol','*','Widths',0.8);

xlabel("Age (Years)",'FontSize',20,"FontName","Times",'FontWeight','bold')
ylabel("Delay (Seconds)",'FontSize',20,"FontName","Times",'FontWeight','bold')

picturewidth = 20; 
hw_ratio = 0.65;
set(findall(hfig,'-property','Box'),'Box','off')
pos = get(hfig,'Position');
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
set(gca,'fontsize',14)
if type == "HR" && nm == 690
    title("HR 690 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 0.5])
elseif type == "HR" && nm == 830
    title("HR 830 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 0.5])
elseif type == "RR" && nm == 690
    title("RR 690 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 2.5])
else
    title("RR 830 nm",'FontSize',18,"FontName","Times",'FontWeight','bold')
    ylim([0 2.5])
end