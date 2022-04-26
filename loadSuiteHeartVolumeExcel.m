function [volume_data] = loadSuiteHeartVolumeExcel(dirPath, excelFile)


if nargin < 1 % Select data
    [excelFile,dirPath] = uigetfile('*.*','Select Suite Heart Excel File');
end

%load volume data from excel
[~,TXT]=xlsread([dirPath, excelFile],2);
i=strfind(TXT,'LAx LV Chamber Volumes'); j=strfind(TXT,'LAx LV Endo Contour Volumes');
i=find(not(cellfun('isempty',i))); j=find(not(cellfun('isempty',j)));
volume=str2double(TXT(i+2:j-2,3));

t_volume=str2double(TXT(i+2:j-2,2))/1000; %in seconds
t_volume=t_volume-t_volume(1); %make sure time vector starts at zero

volume_data=[volume, t_volume];