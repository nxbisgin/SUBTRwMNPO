function [ixExclusionSensors,hc_position,ix_hc_closeSara]= get_exclusionAreaBasedHC(fi,R,channel_labels,rad,user_hc_position)

% get the head coil
% fi : path of the files .ds
% rad: radio
% R  : sensors positions
% channel labels:
% output:
% ixExclusionSensors: index of the SARA sensors belong to exclusion area 
% hc_position: head coil position
% ix_hc_closeSara :   index of the SARA sensor close to head coil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%        Diana Escolana Vargas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    rad=10;
else
    rad=rad;
end

if nargin>=5
    if iscellstr(user_hc_position) || ischar(user_hc_position)
       ix=FindChannels(channel_labels,{user_hc_position});
       user_hc_position=R(ix,:);
    else
       user_hc_position=user_hc_position;
    end
end
[~,~,ext]=fileparts(fi);
if strncmp('.ds',ext,3) && nargin<5
    %%% First get data info with sensor structures
    ctf = readCTFds (fi);
    % extract the hc position
    hc_position= [ctf.hc.dewar(1,4); ctf.hc.dewar(2,4); ctf.hc.dewar(3,4)];
else
    hc_position=user_hc_position;
end
    

% Calculate wich channel is closest to the hc
distance_hc =zeros(1,length(R));

for i=1:length(R)
    distance_hc(1,i) = sqrt((R(i,1)-hc_position(1))^2+(R(i,2)-hc_position(2))^2+(R(i,3)-hc_position(3))^2);
end
[~,ix_hc_closeSara]=min(distance_hc);


for i=1:length(R)
    distance_hc_Sara(1,i) = sqrt((R(i,1)-R(ix_hc_closeSara,1))^2+(R(i,2)-R(ix_hc_closeSara,2))^2+(R(i,3)-R(ix_hc_closeSara,3))^2);
end
% get data around rad , 10cm defauld 
[~,ixExclusionSensors]=find(distance_hc_Sara<=rad);
   
str=sprintf ('Channel closest to the hc : %s',channel_labels{ix_hc_closeSara,:});
disp(str);

return
