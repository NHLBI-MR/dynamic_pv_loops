function data = processPVLoops(prime_data, volume_data,ecg_channel, ecgThresh)

if nargin < 3
    ecg_channel=2;
    ecgThresh = 0.7; %SOMETIMES THIS PARAMETER MUST BE TWEAKED for R peak detection
end

%unpack prime signals
[lvp, aop]=primeAoLV(prime_data(:,end-2), prime_data(:,end-1));
ecg=prime_data(:,ecg_channel);
t_pressure=prime_data(:,end); %time pressure and ecg volume signals
%unpack volume signals
volume = volume_data(:,1);
t_volume = volume_data(:,2); %time vector volume signal

%detect MR gradient activity
gradient_inds = findScanTime(prime_data(:,end-5:end-3)); %gradient channels

%crop, sync and down sample pressures base on detected gradient activity and t_volume
downsampligfactor=4;
[t_pressure, lvp, aop, ecg, t_volume, volume] = reformatSignals(gradient_inds, t_pressure, lvp, aop, ecg, t_volume,volume, downsampligfactor);

%Detect end diastole and end systolce in pressure curve
[ED_inds, ES_inds, detectionsED, detectionsES] = indexEDESpressure(t_pressure, lvp, ecg, ecgThresh);

%find ED and ES in volume data
[ED_volume_inds, ES_volume_inds, volume_detectionsED, volume_detectionsES, volume] = indexEDESvolume(t_volume, volume, t_pressure, ED_inds);

%Plot physiological signals and ED and ES detections
plotSignalsAndDetections(t_pressure, lvp, aop, ecg, t_volume, volume,  detectionsED, detectionsES, volume_detectionsED, volume_detectionsES);

%synchronize pressure and volume data for each individual heart beat
[data] = syncPV_EDES(t_volume, volume, t_pressure, lvp, aop, ED_volume_inds, ES_volume_inds, ED_inds, ES_inds);

%store more data in struct
data(1).Pressure=lvp;
data(1).TPressure=t_pressure;
data(1).Volume=volume;
data(1).TVolume=t_volume;
data(1).ECG=ecg;
data(1).ED_inds=ED_inds;
data(1).ES_inds=ES_inds;
data(1).ED_volume_inds=ED_volume_inds;
data(1).ES_volume_inds=ES_volume_inds;

% %Plot individual PV Loops
data = plotIndividualBeatsAndLoops(data,[]);

% calc ESPVR and EDPVR
data = calcPVRS(data);

%Plot all PV loops
plotAllPVLoops(data);

disp('Finished. Results are stored in data struct.')
end

%----SUPPORT FUNCTIONS----

function [lvp, aop]=primeAoLV(ibp0_data, ibp1_data)
% Determine which signal is LV and aortic pressure

ibp0_scaled = ibp0_data;
ibp1_scaled = ibp1_data;

if (max(ibp0_scaled)-min(ibp0_scaled)) > (max(ibp1_scaled)-min(ibp1_scaled))
    lvp = ibp0_scaled;
    aop = ibp1_scaled;
else
    aop = ibp0_scaled;
    lvp = ibp1_scaled;
end

end

%------------
function [gradient_inds] = findScanTime(gradients)
%gradients: matrix with up to 3 gradient signals
% gradient_inds: detected indeces where MRI scanning was performed

scanIndex=zeros(size(gradients));


for k=1:size(gradients,2) %loop over three gradient recordings
    g=gradients(:,k)-mean(gradients(:,k)); %subtract mean gradient to center singal around zero
    thresh=2*max(abs(g(1:1000))); %threshold is double the max value during first second of recording (no gradients on)
    ind=find(abs(g)>thresh);
    scanIndex(ind,k)=1;
end

scanIndex=sum(scanIndex')>0;

gradient_inds=[];
count=1;
nbrOfPoints=5;
for i =nbrOfPoints+1:round(length(scanIndex)*0.95)
    if  scanIndex(i) == 1 &&(sum(scanIndex(i-nbrOfPoints:i-1)) == 0)
        gradient_inds(count,1)=i;
    end
    if  scanIndex(i) == 1 && (sum(scanIndex(i+1:i+nbrOfPoints)) == 0)
        gradient_inds(count,2)=i;
        count=count+1;
    end
end

if any(gradient_inds(:,2)-gradient_inds(:,1))<10000 %remove to short intervals
    row=find((gradient_inds(:,2)-gradient_inds(:,1))<10000);
    gradient_inds(row,:)=[];
end
end


%------------
function [t_pressure, lvp, aop, ecg, t_volume,volume] = reformatSignals(gradient_inds, t_pressure, lvp, aop, ecg, t_volume,volume, downsampligfactor)

%crop signals to detected gradient activity
t_pressure = t_pressure(gradient_inds(1,1):gradient_inds(1,2))-t_pressure(gradient_inds(1,1)); %time starting at zero
lvp=lvp(gradient_inds(1,1):gradient_inds(1,2));
aop=aop((gradient_inds(1,1):gradient_inds(1,2)));
ecg=ecg(gradient_inds(1,1):gradient_inds(1,2));

%sync time vector for pressue data
totalScanTime=t_volume(end);
startMRI_t=(t_pressure(end)-totalScanTime);
start_ind=find(abs(t_pressure-startMRI_t)==min(abs(t_pressure-startMRI_t)));
t_pressure=linspace(0,t_volume(end), length(t_pressure(start_ind:end)))';
lvp=lvp(start_ind:end); aop=aop(start_ind:end); ecg=ecg(start_ind:end);

%down sample pressure data with a factor
lvp=interp1(t_pressure, lvp, linspace(t_pressure(1), t_pressure(end), length(lvp)/downsampligfactor));
aop=interp1(t_pressure, aop, linspace(t_pressure(1), t_pressure(end), length(aop)/downsampligfactor));
ecg=interp1(t_pressure, ecg, linspace(t_pressure(1), t_pressure(end), length(ecg)/downsampligfactor));
t_pressure=linspace(t_pressure(1),t_pressure(end), length(lvp))';

%interpolate volume to same length as the pressure
volume=interp1(t_volume, volume,linspace(t_volume(1), t_volume(end), length(lvp)))'; %interpolate to length of pressure data
t_volume=linspace(t_volume(1),t_volume(end), length(volume))';

end

%------------
function [ED_inds, ES_inds, detectionsED, detectionsES] = indexEDESpressure(t_pressure, lvp, ecg, ecgThresh)

%normalize ecg from 0 - 1
minECG = min(ecg);
ecg_normalized = ecg -  minECG;
maxECG = max(ecg_normalized);
ecg_normalized = ecg_normalized./maxECG;
[~,locs] = findpeaks(ecg_normalized,t_pressure,'MinPeakHeight',ecgThresh,'MinPeakDistance',0.150); %detect r peaks in ecg
detectionsED=locs';


dt=t_pressure(2)-t_pressure(1);
%loop over detected R peaks to find end diastolic and end systolic indicies
for i=1:length(detectionsED)
    ED_inds(i)=find(t_pressure==detectionsED(i));
    
    if i>1 %end systole: search for dP/dt min and move back 40 ms
        buffer=round(length(ED_inds(i-1):ED_inds(i))/4);
        dP=gradient(lvp(ED_inds(i-1)+buffer:ED_inds(i)), dt);
        IVRT=find(dP==min(dP)); %find IVRT
        ES_inds(i-1)=IVRT+ED_inds(i-1)+buffer-round(0.04/mean(diff(t_pressure)));%index of ES = IVRT-40 ms base on IVRT curation typically 80 ms
    end
end

detectionsES= t_pressure(ES_inds)';
end

%------------
function [ED_volume_inds, ES_volume_inds, volume_detectionsED, volume_detectionsES, volume] = indexEDESvolume(t_volume, volume, t_pressure, detections)


manualCorr=1;

ED_ind_pressure=detections;
[~, ED_ind_pressure] = min(abs(t_volume-t_pressure(detections(1))));

windowED=100;

%circshift to align first ED in volume to the pressure detection
if (ED_ind_pressure(1)-windowED)<1
    windowED=100+(ED_ind_pressure(1)-windowED)-1;
end
[~,firstED]=max(volume(ED_ind_pressure(1)-windowED:ED_ind_pressure(1)+windowED));
firstED=ED_ind_pressure(1)-windowED+firstED-1;
shiftedED=firstED-ED_ind_pressure(1);
volume=circshift(volume,-shiftedED);
firstED=ED_ind_pressure(1);

ED_volume_inds(1)=firstED;

%find ED in pressure time stamp in volume signal
for i=2:length(detections)
    [~, ED_ind_pressure(i)] = min(abs(t_volume-t_pressure(detections(i))));
end


%detect ED as maximum volume around the detected ED in the pressure signal
windowED=70;
[volume_detectionsED, ED_volume_inds] = findEDIndex(t_volume, volume, ED_ind_pressure, windowED);

%re-iterate detection with more narrwo window
windowED=20;
[volume_detectionsED , ED_volume_inds] =findEDIndex(t_volume, volume, ED_volume_inds, windowED);



%show detections and allow manual corrections
figure(222); clf; hold on;
plot(t_volume,volume);
legend('LV volume (ml)','AutoUpdate','off')
xlabel('Time (s)');
ylabel('LV volume (ml)')

%plot detection of each beat
for i=1:length(volume_detectionsED)
    plot([volume_detectionsED(i), volume_detectionsED(i)], [min(volume(:)),max(volume(:))], 'k-', 'LineWidth', 1);
end




display('MOVE VOLUME ED DETECTIONS')
fix = 0;

while(fix<1)
    d = input('Move any detections?   (1/0)           ');
    if(d == 0)
        fix = 1;
        break;
    else
        
        %Manual corrections of ED detections
        disp('Click on detection to move.');
        [x,y]=ginput(1);
        [~, moveInd] = min(abs(volume_detectionsED-x));
        
        %UPDATE PLOT
        figure(222); hold off;
        plot(t_volume,volume);
        hold on;
        for i=[1:moveInd-1, moveInd+1:length(volume_detectionsED)]
            plot([volume_detectionsED(i), volume_detectionsED(i)], [min(volume(:)),max(volume(:))], 'k-', 'LineWidth', 1);
        end
        
        disp('Click on new location.');
        [x,y]=ginput(1);
        [~, newInd] = min(abs(t_volume-x));
        windowED=10;
        [newT , newInd] = findEDIndex(t_volume, volume, [ED_volume_inds(1), newInd], windowED);
        ED_volume_inds(moveInd)=newInd(2);
        volume_detectionsED(moveInd)=newT(2);
        
        
        
        %UPDATE PLOT
        figure(222); clf; hold on;
        plot(t_volume,volume);
        legend('LV volume (ml)','AutoUpdate','off')
        xlabel('Time (s)');
        ylabel('LV volume (ml)')
        
        %plot detection of each beat
        for i=1:length(volume_detectionsED)
            plot([volume_detectionsED(i), volume_detectionsED(i)], [min(volume(:)),max(volume(:))], 'k-', 'LineWidth', 1);
        end
    end
end


volume_detectionsED = unique(sort([t_volume(ED_volume_inds)']));
ED_volume_inds= unique(sort(ED_volume_inds));

%detect ES as minimum pressure in first part of each detected cardiac cycle
for i = 2:length(ED_volume_inds)
    searchRange=ED_volume_inds(i-1):ED_volume_inds(i);
    searchRange=searchRange(2:round(length(searchRange)*2/3));
    [~,volind]=min(volume(searchRange));
    ES_volume_inds(i-1)=volind(1)+ED_volume_inds(i-1)-1;
end

volume_detectionsES = t_volume(ES_volume_inds)';

close(222);
end


%------------
function [volume_detectionsED, ED_volume_inds] =findEDIndex(t_volume, volume, ED_ind_pressure, windowED)

ED_volume_inds(1)=ED_ind_pressure(1);

%detect ED as maximum volume around the detected ED in the pressure signal
for i=2:length(ED_ind_pressure)
    
    %define search range for ED in volume
    if ED_ind_pressure(i)-windowED<1
        searchRange =1:ED_ind_pressure(i)+windowED;
    elseif ED_ind_pressure(i)+windowED>length(volume)
        searchRange =ED_ind_pressure(i)-windowED:length(volume);
    else
        searchRange = ED_ind_pressure(i)-windowED:ED_ind_pressure(i)+windowED;
    end
    
    %don't allow a search region overlapping with the previous detection
    if any(searchRange==ED_ind_pressure(i-1))
        ii=find(searchRange==ED_ind_pressure(i-1));
        if ((ii+windowED)>length(searchRange))
            searchRange = searchRange(length(searchRange)-20):searchRange(end);
        else
            searchRange = searchRange(ii+windowED):searchRange(end);
        end
    end
    
    [~,volind]=max(volume(searchRange));
    ED_volume_inds(i)=searchRange(1)+volind(1)-1;
end


volume_detectionsED = unique(sort([t_volume(ED_volume_inds)']));

end

%------------
function plotSignalsAndDetections(t_pressure, lvp, aop, ecg, t_volume, volume,  detections, detectionsES, volume_detectionsED, volume_detectionsES)


%Plot physiological signals
figure(8); clf;

%ECG
subplot(3,1,1); hold on;
plot(t_pressure,ecg);
legend('ECG (ml)','AutoUpdate','off')
xlabel('Time (s)');
ylabel('ECG')

%LV and Ao pressures
subplot(3,1,2); hold on;
plot(t_pressure,lvp,t_pressure,aop);
legend('LV pressure (mmHg)', 'Aortic pressure (mmHg)','AutoUpdate','off')
xlabel('Time (s)');
ylabel('Pressure (mmMg)')

%LV volumes
subplot(3,1,3); hold on;
plot(t_volume,volume);
legend('LV volume (ml)','AutoUpdate','off')
xlabel('Time (s)');
ylabel('LV volume (ml)')



%plot detection of each beat
for i=1:length(detections)
    
    subplot(3,1,1); hold on;
    plot([detections(i),detections(i)], [min(ecg(:)),max(ecg(:))], 'k:', 'LineWidth', 0.1);
    plot([detections(i)], ecg(detections(i)==t_pressure), 'ko');
    
    subplot(3,1,2); hold on;
    plot([detections(i),detections(i)], [min([lvp(:); aop(:)]),max([lvp(:); aop(:)])], 'k-', 'LineWidth', 1);
    
    subplot(3,1,3); hold on;
    plot([volume_detectionsED(i), volume_detectionsED(i)], [min(volume(:)),max(volume(:))], 'k-', 'LineWidth', 1);
    
    if i>1
        subplot(3,1,1); hold on;
        plot([detectionsES(i-1), detectionsES(i-1)],  [min(ecg(:)),max(ecg(:))], 'k:', 'LineWidth', 1)
        
        subplot(3,1,2); hold on;
        plot([detectionsES(i-1), detectionsES(i-1)],  [min([lvp(:); aop(:)]),max([lvp(:); aop(:)])], 'k:', 'LineWidth', 1); %plot detection of each beat
        
        subplot(3,1,3); hold on;
        
        plot([volume_detectionsES(i-1), volume_detectionsES(i-1)], [min(volume(:)),max(volume(:))], 'b:', 'LineWidth', 1);
        
    end
end
end

%------------
function [data, ES] = syncPV_EDES(t_volume, volume, t_pressure, lvp, aop,  ED_volume_inds, ES_volume_inds, ED_inds, ES_inds)


figure(1000); clf;
for i=2:length(ED_volume_inds)
    
    %extract pressure curve and find dicrotic notch
    beatInds=ED_inds(i-1):ED_inds(i);
    ao=aop(beatInds);
    p=lvp(beatInds);
    
    EST_pressure=find(p==lvp(ES_inds(i-1)));
    
    
    %extract volume curve and find end systole
    beatIndsV=ED_volume_inds(i-1):ED_volume_inds(i);
    v=volume(beatIndsV);
    t_v=linspace(0,1,length(v));
    
    EST=find(v==volume(ES_volume_inds(i-1)));
    %     [~,EST]=min(v);
    
    
    %sync volume curve to pressure curve
    v_sys=v(1:EST);
    t_sys=t_v(1:EST);
    v_sys=interp1(t_sys, v_sys ,linspace(t_sys(1), t_sys(end), EST_pressure)); %interpolate to length of pressure data
    
    v_dias=v(EST:end);
    t_dias=t_v(EST:end);
    v_dias=interp1(t_dias, v_dias ,linspace(t_dias(1), t_dias(end), length(p)-EST_pressure+1));
    v=[v_sys, v_dias(2:end)];
    
    tRR=linspace(0, t_pressure(ED_inds(i))-t_pressure(ED_inds(i-1)), length(p));
    
    %find upper corner in PV loop
    dt= t_pressure(2)-t_pressure(1);
    ES = findCorner(p,v',dt);
    
    %save ESP in multiple pv loop arrays as well
    data(i-1).lvp=p(:)';
    data(i-1).aop=ao(:)';
    data(i-1).lvv=v(:)';
    data(i-1).t=tRR(:)';
    data(i-1).ES=ES';
    data(i-1).HR=round(60./tRR(end));
    data(i-1).EDV=v(1);
    data(i-1).ESV=v(ES);
    data(i-1).ESP=p(ES);
end


end

%------------
function [ESP] = findCorner(p,v,dt)

v=v(:); p=p(:);
%cornerPoint detection according to appendix in Lankhaar et al 2009
%Modeling the Instantaneous Pressure–Volume Relation of the Left Ventricle: A Comparison of Six Models

dP=gradient(p,dt);
dPdtmin_ind=find(dP==min(dP));
dPdtmin_ind=dPdtmin_ind(1);

minP=min(p);
maxP=max(p);

h=(maxP-minP); %height
Pcenter=minP+ h/2; %center point

minV=min(v(1:dPdtmin_ind));
maxV=v(1);

w=(maxV-minV); %widht
Vcenter=minV+w/2; %center point

searchInds = find((v'<Vcenter)+(p'>Pcenter) ==2);
d=sqrt(((v(searchInds)-Vcenter)/w).^2 + ((p(searchInds)- Pcenter)/h).^2);
[~,ii] = max(d);
ESP=searchInds(ii);

end


%------------
function data = plotIndividualBeatsAndLoops(data, firstBeat, nbrOfBeats)

if nargin < 3
    nbrOfBeats=10;
end
if nargin<4
    fignbr=110;
end

figure(fignbr); clf;
subplot(ceil(sqrt(length(data)))+2, ceil(sqrt(length(data))), [1:ceil(sqrt(length(data)))] ); hold on;
plot(data(1).TVolume, data(1).Volume)
plot(data(1).TVolume(data(1).ED_volume_inds), data(1).Volume(data(1).ED_volume_inds),'ro')
plot(data(1).TVolume(data(1).ES_volume_inds), data(1).Volume(data(1).ES_volume_inds),'bo')
subplot(ceil(sqrt(length(data)))+2, ceil(sqrt(length(data))), [ceil(sqrt(length(data)))+1:ceil(sqrt(length(data)))*2] ); hold on;
plot(data(1).TPressure, data(1).Pressure)
plot(data(1).TPressure(data(1).ED_inds), data(1).Pressure(data(1).ED_inds),'ro')
plot(data(1).TPressure(data(1).ES_inds), data(1).Pressure(data(1).ES_inds),'bo')

for i=1:length(data)
    subplot(ceil(sqrt(length(data)))+2, ceil(sqrt(length(data))), ceil(sqrt(length(data)))*2+i); hold on;
    p=data(i).lvp;
    v=data(i).lvv;
    ES = data(i).ES;
    
    plot(v,p ,'k') %PV loop
    plot(v(1), p(1),'ro') %End diastole
    plot(v(ES), p(ES),'bo') %End systole
    
    title(num2str(i))
    xlabel('LV volume (ml)');
    ylabel('LV pressure (mmHg)')
    xlim([min([data.lvv])-10, max([data.lvv])+10])
    ylim([min([data.lvp])-10, max([data.lvp])+10])
    axis square
end


if isempty(firstBeat)
    disp('SELECT PV LOOPS TO USE');
    input_beat = input('First of 10 PV loops to use   (0 if all loops.)           ');
else
    input_beat=firstBeat;
end

if input_beat~=0
    d=[1:input_beat-1, input_beat+nbrOfBeats:length(data)]; %data points to remove
else
    d=0;
end

if d~=0
    data(d)=[]; %remove non-selected beats
    
    figure(fignbr); clf;
    for i=1:length(data)
        subplot(ceil(sqrt(length(data))), ceil(sqrt(length(data))), i); hold on;
        p=data(i).lvp;
        v=data(i).lvv;
        ES = data(i).ES;
        
        %     plot([v, v(1)],[p, p(1)] ,'k') %PV loop
        plot(v,p ,'k') %PV loop
        plot(v(1), p(1),'ro') %End diastole
        plot(v(ES), p(ES),'bo') %End systole
        
        title(num2str(i))
        xlabel('LV volume (ml)');
        ylabel('LV pressure (mmHg)')
        
        xlim([min([data.lvv])-10, max([data.lvv])+10])
        ylim([min([data.lvp])-10, max([data.lvp])+10])
        axis square
    end
    %     legend('PV loop', 'End diastole','End systole', 'Location', 'NorthEastOutside','FontSize', 12)
    
end

if isempty(firstBeat)
    
    display('ADJUST ESPVR POINT DETECTIONS?')
    d = input('Loops to adjust   (0 if none.)           ');
else
    d=0;
end

if d~=0
    for loop=d
        %Manual corrections of ED detections
        disp('Click on new corner point.');
        [x,y]=ginput(1);
        [~, ES] = min(abs(data(loop).lvv-x) + abs(data(loop).lvp-y));
        plot(data(loop).lvv(ES),data(loop).lvp(ES),'bo', 'MarkerSize', 16)
        data(loop).ES=ES;
    end
    
    
    figure(fignbr); clf;
    for i=1:length(data)
        subplot(ceil(sqrt(length(data))), ceil(sqrt(length(data))), i); hold on;
        p=data(i).lvp;
        v=data(i).lvv;
        ES = data(i).ES;
        
        plot(v,p ,'k') %PV loop
        plot(v(1), p(1),'ro') %End diastole
        plot(v(ES), p(ES),'bo') %End systole
        
        title(num2str(i))
        xlabel('LV volume (ml)');
        ylabel('LV pressure (mmHg)')
        xlim([min([data.lvv])-10, max([data.lvv])+10])
        ylim([min([data.lvp])-10, max([data.lvp])+10])
        axis square
    end
    
    
end

data(1).InputBeat=input_beat;

end

%------------
function [data] = calcPVRS(data)

for i=1:length([data.ES])
    p=data(i).lvp;
    v=data(i).lvv;
    ES = data(i).ES;
    
    volES(i)=v(ES);
    lvpES(i)=p(ES);
    
    volED(i)=v(1);
    lvpED(i)=p(1);
end

[volES, inds]=sort(volES);
lvpES=lvpES(inds);

[volED, inds]=sort(volED);
lvpED=lvpED(inds);

%linear regerssion for ESPVR (ESP = k*ESV + m)
ESPVR_f = polyfit(volES,lvpES,1); %first degree polynomial
V0= -ESPVR_f(2)/ESPVR_f(1); %analytically derive V0
ESPVR = polyval(ESPVR_f,sort([volES]));

%linear regerssion for ESPVR (EDP = k*EDV + m)
EDPVR_f = polyfit(volED,lvpED,1); %first degree polynomial
EDPVR = polyval(EDPVR_f,sort([volED]));

data(1).V0=V0;
data(1).Contractility=ESPVR_f(1);
data(1).ESPVR=ESPVR;
data(1).Compliance=1./EDPVR_f(1);
data(1).EDPVR=EDPVR;

data(1).ESPVR_f=[ESPVR_f(1), ESPVR_f(2)];
data(1).EDPVR_f=[EDPVR_f(1)  EDPVR_f(2)];
end
%------------
function plotAllPVLoops(data)

for i = 1:length(data)
    
    p=data(i).lvp;
    v=data(i).lvv;
    ES = data(i).ES;
    
    volES(i)=v(ES);
    lvpES(i)=p(ES);
    
    volED(i)=v(1);
    lvpED(i)=p(1);
    
end
%Plot all PV loops
figure; clf; hold on;

plot([data(1).lvv, data(1).lvv(1)],[data(1).lvp, data(1).lvp(1)] ,'k') %PV loop
plot(sort([volES]), data(1).ESPVR,'b', 'LineWidth', 2)
plot(sort([volED]), data(1).EDPVR,'r', 'LineWidth', 2)
legend('PV loops', 'ESPVR','EDPVR', 'Location', 'NorthEastOutside','FontSize', 14, 'AutoUpdate', 'Off')
xlabel('LV volume (ml)');
ylabel('LV pressure (mmHg)')

plot(volED, lvpED, 'ro', 'LineWidth', 2); %ES
plot(volES, lvpES, 'bo', 'LineWidth', 2); %ED

for i = 1:length(data)
    p=data(i).lvp;
    v=data(i).lvv;
    plot([v, v(1)],[p, p(1)] ,'k:') %PV loop
end
end