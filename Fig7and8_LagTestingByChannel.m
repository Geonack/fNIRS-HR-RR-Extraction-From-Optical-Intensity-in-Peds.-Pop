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
if type == "HR";
    load("AbsoluteValueLagsHR.mat")
else
    load("AbsoluteValueLagsRR.mat")
end

plotlags = meanlags';
if nm == 690
    ranges = 1:46;
else
    ranges= 47:92;
end
newMatrix = mean(plotlags(:,ranges));
newMatrix2 = mean(plotlags(:,1:92));


X1 = 1:47;
%% Final Figure Histogram Plot by channel
clc;

figure1 = figure();
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot1 = boxchart(plotlags(:,ranges),"Notch","on");
hold on
bar(X1,[newMatrix mean(newMatrix)],"FaceColor",[0.501960784313725 0.501960784313725 0.501960784313725]);
alpha(.5);
hold off
ylabel('Delay (Seconds)','FontWeight','bold','FontName','Times New Roman','FontSize',25);

xlabel('Channels','FontWeight','bold','FontName','Times New Roman','FontSize',25);
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

%% Channel Delay Mapping on image SC
channelMap = [0 0 1 2 0 5 6 0 0;0 0 3 4 0 7 8 0 0; 9 10 13 14 0 17 18 21 22;11 12 15 16 0 19 20 23 24;...
   25 26 28 29 0 32 33 36 37; 0 27 30 31 0 34 35 38 0; 0 0 29 40 0 43 44 0 0;0 0 41 42 0 45 46 0 0 ];
dispMap = NaN([8 9]);
for i = 1:46
    loc = find(channelMap==i);
    dispMap(loc) = newMatrix(i);
end
temp = reshape(dispMap,[1 72]);
temp = normalize(temp);%,"range");
temp = reshape(temp,[8 9]);

figure;
h = heatmap(temp,"FontSize",15,"FontName","Times")
h.Colormap = turbo(256)
clim([-1.9962 3.7558])
h.MissingDataColor = [1 1 1]
grid off
max(temp,[],"all")
min(temp,[],"all")

mean(dispMap,'all',"omitnan")
std(dispMap,[],"all","omitnan")
    
%% New Analysis
clearvars -except plotlags newMatrix ranges nm type
formatSpec = 'Channel%d';
for A1 = 1:46
    channames(A1) = string(sprintf(formatSpec,A1));
end
lags690 = array2table(plotlags(:,ranges),'VariableNames',channames);
lags690 = addvars(lags690,[1:107]','Before','Channel1','NewVariableNames',"Subjects");
Meas = table([1:46]',VariableNames="Channels");
rm = fitrm(lags690,'Channel1-Channel46~1','WithinDesign',Meas)
%rm = fitrm(lags690,'Channel1-Channel46~Subjects','WithinDesign',Meas)

ranovatbl = ranova(rm)
spheric = mauchly(rm)
epsilons = epsilon (rm)
comp = multcompare(rm,"Channels",'ComparisonType','lsd');
comp_matrix = NaN(46);
for i = 1:length(comp.Channels_1)
    comp_matrix(comp{i,1},comp{i,2})=comp{i,5};
end

for i = 1:46
    for j = 1:46
      [~,~,~,stats] = ttest(lags690{:,i+1},lags690{:,j+1});
      tstatistic(i,j) = stats.tstat;
    end
end
tstatistic = abs(tstatistic);
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(comp_matrix,.05,'dep','yes');

[x1,y1] = find(h);

figure1 = figure();
axes1 = axes('Parent',figure1);
imagesc(tstatistic);
colorbar;
hold on
for i = 1:length(x1)
    rectangle('Position', [x1(i)-0.5,y1(i)-0.5,1,1],'LineWidth',1,...
        'Curvature',[1,1],"EdgeColor","Red")  % Highlight significant ones (after FDR) with curved rectangles (basically circles)
end
hold off
ylabel('Channels','FontWeight','bold','FontName','Times New Roman','FontSize',25);
xlabel('Channels','FontWeight','bold','FontName','Times New Roman','FontSize',25);
if type == "HR" && nm == 690
    title("HR 690 nm")
elseif type == "HR" && nm == 830
    title("HR 830 nm") 
elseif type == "RR" && nm == 690
    title("RR 690 nm") 
else
    title("RR 830 nm") 
end
box(axes1,'on');
hold(axes1,'off');
xlim([0.5 46.5]);
ylim([0.5 46.5]);
xticks([5:5:45]);
yticks([5:5:45]);
if type == "HR"
    clim([0 14])
else
    clim([0 8.5])
end
set(axes1,'FontName','Times New Roman','FontSize',15);
set(findall(figure1,'-property','Box'),'Box','off') 
set(findall(figure1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figure1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis square
length(x1)