clear;
close all;
load T1fit;
load T2fit;


mT1est = T1est;

T1s = 1000*bsxfun(@times, T1est, mask);
% 
% h=figure 
% 
% imshow(T1s, [0, 2500]);
% colormap('parula'), colorbar;
% title('T1 map (ms)', 'FontSize', 22); 
% %faxis
%%
x=figure

hist(mT1est(mT1est>0.05 & mT1est < 2.5)*1000, 400); 
%faxis

title('T1 Histogram', 'FontSize', 22)
xlabel('T1 (ms)')
ylabel('Voxel Count')
axis([00 2500 0 1250])
set(gca,'FontSize', 20)

%%


mT2est = T2est;
T2s = 1000*bsxfun(@times, mT2est, T2mask);
% y=figure
% 
% imshow(T2s, [0, 160]); colormap('parula'), colorbar;
% title('T2 map (ms)', 'FontSize', 22); %faxis

%%

z=figure

mT2vals = mT2est(mT2est~=0);
hist(mT2vals(mT2vals>0.006 & mT2vals<.25)*1000, 400); %faxis
title('T2 Histogram', 'FontSize', 22)
xlabel('T2 (ms)')
ylabel('Voxel Count')
axis([0 250 0 1250])
set(gca,'FontSize', 20)

% figure
% scat_mask = T2est < .300 & T1est < 3.5;
% plot(reshape(T2est(scat_mask),1,[]),reshape(T1est(scat_mask),1,[]), 'o')
%%

% figure
% T2vals = mT2est(mT2est~=0);
% figure; hist(T2vals(T2vals<.3)*1000, 100); faxis
% 
% figure;
% scat_mask = T2est > 0.01 & T2est < .300 & T1est > 0.01 & T1est < 3.5;
% hist3([reshape(1000*T1est(scat_mask),[],1), reshape(1000*T2est(scat_mask),[],1)],'Nbins',[100,100], 'CDataMode','auto','FaceColor','interp')