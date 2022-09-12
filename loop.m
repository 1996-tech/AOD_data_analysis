clc
clear all
close all

load('workspace_dec21-2nd16.mat')

for i=1:16
    GridID=zeros(426277,1);
    GridID(:)=[1:426277];
    latitude=[m_coords(:,1)];
    longitude=[m_coords(:,2)];
    timestamp=[zeros(426277,1)];
    timestamp(:,1)=[v_time(i)];
    Date=datetime(timestamp,'ConvertFrom','dateNum','Format','yyyy-MM-dd');
    AOD1=[transpose(a_data(i,:,1))];
    AOD2= [transpose(a_data(i,:,2))];
    Month=month(Date);
    Week=week(Date);
    Jan2020=table(GridID,latitude,longitude,timestamp,Date,AOD1,AOD2,Month,Week);
    writetable(Jan2020,['D' num2str(i) '_' num2str(Month(1)) '_' '2021.csv'])
end

