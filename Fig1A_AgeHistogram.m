% Figure 1 Age Histogram code
clc;
clear;
close all;

% Choose Directory for Subject Data
cd %CD subject folder directory

%Load Subject Folder
allsubject= dir('*.snirf');
subjects = 1:107; 
peakcutoff = .9;
%% Load Raw Data
parfor subjnumber=subjects     
    %disp("Analyzing subject number: "+subjnumber)
    if subjnumber<10
        allsubjectidx=subjnumber;
    elseif subjnumber>99
        allsubjectidx=subjnumber-90;
    else
        allsubjectidx=subjnumber+8;
    end

    subjname=allsubject(allsubjectidx).name;
    raw=SnirfLoad(subjname);
    
    %% Add age to data

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
end 

%% Create figure
hfig = figure;
axes1 = axes('Parent',hfig);
hold(axes1,'on');
histogram(agedata,'Parent',axes1,...
    'FaceColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'BinMethod','auto');
ylabel('Number of Subjects','FontWeight','bold','FontName','Times New Roman','FontSize',20);
xlabel('Age','FontWeight','bold','FontName','Times New Roman','FontSize',20);
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontName','Times New Roman','MinorGridLineStyle','-','YMinorGrid','on');
picturewidth = 20; 
hw_ratio = 0.65;
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
axis("square")