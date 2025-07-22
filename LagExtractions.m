clc;
clear;
close all;

%Used for figures 6,7, and 8
load("ChanByChanHRCorr.mat"); %Load data saved from figure 4 script

%Uncomment whether HR or RR lags should be extracted
type = "HR"
%type = "RR"

if type == "HR"
    HRdata = allDataHR(:,:,177);
else
    HRdata = allDataRR(:,:,174);
end

cd %CD subject folder directory

clearvars -except RRdata HRdata allsubject type

for ii = 1:107
    HRanalyzisdata = HRdata(HRdata(:,ii)~=0,ii);
    outHR = isoutlier(HRanalyzisdata,"mean");
    HRanalyzisdata = HRanalyzisdata(~outHR);
    HRmeans(ii) = mean(HRanalyzisdata);
end

for ii = 1:107
    if ii<10
        allsubjectidx=ii;
    elseif ii>99
        allsubjectidx=ii-90;
    else
        allsubjectidx=ii+8;
    end

    subjname=allsubject(allsubjectidx).name;
    raw=SnirfLoad(subjname);
    disp("Analyzing subject number: "+ii)

    if type == "HR"
        HRrawdata = hmrR_BandpassFilt(raw.data,1,2);
    else
        HRrawdata = hmrR_BandpassFilt(raw.data,.2,.5);
    end
    HRmeanband = hmrR_BandpassFilt(raw.data,HRmeans(ii),HRmeans(ii));
    for jj = 1:92
        for kk = 1:92
            [rrcal,lagscal] = xcorr(HRmeanband.dataTimeSeries(:,jj),HRmeanband.dataTimeSeries(:,kk),30,"normalized");
            [pxx,f] = periodogram(rrcal,[],[],50,"onesided");
            [~,locs] = findpeaks(pxx,'SortStr','descend','NPeaks',1);
            [rr,lags] = xcorr(HRmeanband.dataTimeSeries(:,jj),HRmeanband.dataTimeSeries(:,kk),round((1/f(locs))*50),"normalized");
            [Maxi,Indexz] = max(rr);
            lagmat(jj,kk,ii) = abs(lags(Indexz)./50);
        end
    end
    meanlags(:,ii) = mean(lagmat(:,:,ii),2);
end

if type == "HR"
    save("AbsoluteValueLagsHR.mat")
else
    save("AbsoluteValueLagsRR.mat")
end