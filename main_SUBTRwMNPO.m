clear;close all;clc;
addpath(genpath('C:\RESEARCH\NES\SUBTRwMNPOpaper\ctf_software'));
addpath(genpath('C:\RESEARCH\NES\SUBTRwMNPOpaper\MI_functions'));
addpath(genpath('C:\RESEARCH\NES\SUBTRwMNPOpaper\SARA_functions'));
addpath(genpath('C:\RESEARCH\NES\SUBTRwMNPOpaper\MIToolbox'));


%input the data file path
fi='C:\RESEARCH\NES\fBSER\lr28\29w\lr28_AER-STAN_20070402_29w.ds';
[folderpath,fileName] = fileparts(fi);
% read data [samples x channels]
% R: channel positions [#channels x 3]
% sf: sampling frequency
[data,R,sf,channel_names]=getData(fi);
Ns = size(data,1);  % number of samples
Nc = size(data,2);  % number of channels
% find exclusion area, 10 cm around head coil
% ExcCh: Channel numbers for the channels inside exclusion area
% hc_pos: Coordinates of the head coil position that was determined by prior ultrasound
% ch_closestHC: Channel that is closest to the head coil location
radius=10;
[ExcCh,hc_pos,ch_closestHC]= findExclusionArea(fi,R,channel_names,radius);
% HP filter at 0.5 Hz
data = high_pass_filter(0.5,data,sf);

%%%%%%%%%%%%%%%% SUBTRact-mMCG %%%%%%%%%%%%%%%%%%%%
% Find 12 mMCG reference channels based on mRmR on MI. The channels are
% restricted to be outside the exclusion area.
[refCh_mMCG,AmMCG]=find_mMCG_refs(data,sf,channel_names,ExcCh); 
% RefmMCG: Reference signal  from the mMCG reference channels. This is the signal that will be subtracted [chan x samp]
RefmMCG=data(:,refCh_mMCG)';
% Input_mMCG: Input channels for subtraction
% Subtraction will be done on all channels except for mMCG reference channels
% Reference channels after subtraction will become zero
% Excluding them from subtraction will save computational time
Input_mMCG=ones(size(data,2),1);Input_mMCG(refCh_mMCG)=0;Input_mMCG=find(Input_mMCG);
% Raw: The data in input channels for SUBTR [chan x samp]
Raw=data(:,Input_mMCG)';
% apply SUBTRact to remove mMCG
[SUBTR_mMCGtemp]=subtractMethod_optimized(RefmMCG,Raw,sf);
%mMCG reference channels become 0 after SUBTRact, so add them as 0.
SUBTR_mMCG=zeros(size(data))';
SUBTR_mMCG(Input_mMCG,:)=SUBTR_mMCGtemp; % [chan x samp]
   

%%%%%%%%%%%%%%%% SUBTRact-fMCG %%%%%%%%%%%%%%%%%%%%
% Find 18 fMCG reference channels based on mRmR on MI. The channels are
% restricted to be outside the exclusion area and different than the mMCG
% references.
[refCh_fMCG,AfMCG]=find_fMCG_refs(SUBTR_mMCG',sf,channel_names,ExcCh,refCh_mMCG,18);
% ReffMCG: Reference signal  from the fMCG reference channels. This is the signal that will be subtracted [chan x samp]
ReffMCG=SUBTR_mMCG(refCh_fMCG,:);
% Input_fMCG: Input channels for subtraction
% Subtraction will be done on all channels except for mMCG and fMCG reference channels
% Reference channels after subtraction will become zero
% Excluding them from subtraction will save computational time
Input_fMCG=ones(size(data,2),1);Input_fMCG(refCh_fMCG)=0;Input_fMCG(refCh_mMCG)=0;Input_fMCG=find(Input_fMCG);
% Raw_mMCG: The data in input channels for SUBTR [chan x samp]
Raw_mMCG=SUBTR_mMCG(Input_fMCG,:);
% apply SUBTRact for remove fMCG
[SUBTR_fMCGtemp]=subtractMethod_optimized(ReffMCG,Raw_mMCG,sf);
%fMCG reference channels become 0 after SUBTRact, so add them as 0.
SUBTR_fMCG=zeros(size(data))';
SUBTR_fMCG(Input_fMCG,:)=SUBTR_fMCGtemp; % [chan x samp]

% Find mMCG R markers
ulbpm=110;
llbpm=40;
data1.data=data;
data1.Comment=strcat(folderpath,'/',fileName);
data1.sf=sf;
data1.good_channel_labels=channel_names';
[out1]=Umit_Heart_Extract_forMATfilesSARA(data1,ulbpm,llbpm,0,0,'mMCGRmarker');
iRm=out1.indexnew;
clear data1;
    
% Use Orthogonal Projection on Raw data to remove mMCG
OPm = find_orthogonal_projection(data,iRm);
OP_mMCG = data*OPm;
 
% Find fMCG R markers
ulbpm=210;
llbpm=100;
data2.data=OP_mMCG;
data2.Comment=strcat(folderpath,'/',fileName);
data2.sf=sf;
data2.good_channel_labels=channel_names';
[out2]=Umit_Heart_Extract_forMATfilesSARA(data2,ulbpm,llbpm,0,0,'fMCGmarker');
iRf=out2.indexnew;
clear data2;

% Find mMCG only signal 
% Start with raw data 
% Use covariance based Minimum Norm Projection Operator
Lm = time_locked_avg(data,iRm);
Cd = cov(data);
Cm = cov(Lm);
Pm = Cm/Cd;
only_mMCG = data*Pm';

% Find fMCG only signal 
% Start with OP-mMCG 
% Use covariance based Minimum Norm Projection Operator
Lf = time_locked_avg(OP_mMCG,iRf);
Co = cov(OP_mMCG);
Cf = cov(Lf);
Pf = Cf/Co;
only_fMCG = OP_mMCG*Pf';
    
%%%%%%%%%%%%%%% SUBTRwMNPO %%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SUBTRact with Minimum Norm Projection Operator %%%%%%%%%%%%%%%%%
% Input Channels: Channels that SUBTRwMNPO will be performed on
% ExcCh: Exclusion area channels (channels around fetal head)
InCh_mfMCG=ExcCh;
% Input Data: Data to be subtracted 
InData_mfMCG=data(:,InCh_mfMCG)';    %[chan x samp]
% Reference Channels: Channels that will be used as references to SUBTRwMNPO
% ExcCh: Exclusion area channels (channels around fetal head)
RefCh_mfMCG=ExcCh;
% Reference Data: Sum of mMCG only and fMCG only signal from the channels   
RefData_mfMCG=(only_mMCG(:,RefCh_mfMCG)'+only_fMCG(:,RefCh_mfMCG)');
% apply SUBTRwMNPO to remove mMCG and fMCG at the same time
SUBTR_mfMCG=subtractMethod_optimized(RefData_mfMCG,InData_mfMCG,sf);
    
display('Done with SUBTRact');

SUBTR_mMCG=(SUBTR_mMCG)';
SUBTR_fMCG=(SUBTR_fMCG)';
SUBTR_mfMCG=(SUBTR_mfMCG)';

% get samples points to truncate
[stTrunc,endTrunc]=TruncateData(size(SUBTR_fMCG,1),sf,0.4);
% truncate  after SUBTRact
data=data(stTrunc:endTrunc,:);
SUBTR_mMCG=SUBTR_mMCG(stTrunc:endTrunc,:);
SUBTR_fMCG=SUBTR_fMCG(stTrunc:endTrunc,:);
SUBTR_mfMCG=SUBTR_mfMCG(stTrunc:endTrunc,:);
OP_mMCG=OP_mMCG(stTrunc:endTrunc,:);

% get truncated mMCG marker
ulbpm=110;
llbpm=40;
data1.data=data;
data1.Comment=strcat(folderpath,'/',fileName);
data1.sf=sf;
data1.good_channel_labels=channel_names';
[out1]=Umit_Heart_Extract_forMATfilesSARA(data1,ulbpm,llbpm,0,0,'mMCGRmarker2');
iRm2=out1.indexnew;
clear data1;

% get truncated fMCG marker
ulbpm=210;
llbpm=100;
data2.data=OP_mMCG;
data2.Comment=strcat(folderpath,'/',fileName);
data2.sf=sf;
data2.good_channel_labels=channel_names';
[out2]=Umit_Heart_Extract_forMATfilesSARA(data2,ulbpm,llbpm,0,0,'fMCGmarker');
iRf2=out2.indexnew;
clear data2;
    
%filter before average
data = low_pass_filter(40,data,sf);
SUBTR_mMCG = low_pass_filter(40,SUBTR_mMCG,sf);
SUBTR_fMCG = low_pass_filter(40,SUBTR_fMCG,sf);
SUBTR_mfMCG = low_pass_filter(40,SUBTR_mfMCG,sf);
    
display('finding residuals');

% average over mMCG marker after -mMCG
[avg_data_iRm,avg_data_iRm_span]= (time_locked_avg2(data(:,ExcCh), iRm2));
[avg_SUBTR_mMCG_iRm,avg_SUBTR_mMCG_iRm_span ]= (time_locked_avg2(SUBTR_mMCG(:,ExcCh), iRm2));

% average over mMCG marker after -mMCG-fMCG
[avg_SUBTR_fMCG_iRm,avg_SUBTR_fMCG_iRm_span ]= (time_locked_avg2(SUBTR_fMCG(:,ExcCh), iRm2));
[avg_SUBTR_mfMCG_iRm,avg_SUBTR_mfMCG_iRm_span ]= (time_locked_avg2(SUBTR_mfMCG(:,:), iRm2));

% average over fMCG after -mMCG-fMCG
[avg_data_iRf,avg_data_iRf_span]= (time_locked_avg2(data(:,ExcCh), iRf2));
[avg_SUBTR_fMCG_iRf,avg_SUBTR_fMCG_iRf_span ]= (time_locked_avg2(SUBTR_fMCG(:,ExcCh), iRf2));
[avg_SUBTR_mfMCG_iRf,avg_SUBTR_mfMCG_iRf_span ]= (time_locked_avg2(SUBTR_mfMCG(:,:), iRf2));

%%%%%%%%%%%%%%%%%%% find attenuation %%%%%%%%%%%%%%
display('finding attenuation');

[yupper,ylower]=envelope(avg_data_iRm);
mean_env_avg_data_iRm=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_data_iRm(mm)=mean(yupper(mm,:));
end
max_env_avg_data_iRm=max(mean_env_avg_data_iRm); 

[yupper,ylower]=envelope(avg_SUBTR_mMCG_iRm);
mean_env_avg_SUBTR_mMCG_iRm=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_SUBTR_mMCG_iRm(mm)=mean(yupper(mm,:));
end
max_env_avg_SUBTR_mMCG_iRm=max(mean_env_avg_SUBTR_mMCG_iRm); 

[yupper,ylower]=envelope(avg_SUBTR_fMCG_iRm);
mean_env_avg_SUBTR_fMCG_iRm=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_SUBTR_fMCG_iRm(mm)=mean(yupper(mm,:));
end
max_env_avg_SUBTR_fMCG_iRm=max(mean_env_avg_SUBTR_fMCG_iRm); 

[yupper,ylower]=envelope(avg_SUBTR_mfMCG_iRm);
mean_env_avg_SUBTR_mfMCG_iRm=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_SUBTR_mfMCG_iRm(mm)=mean(yupper(mm,:));
end
max_env_avg_SUBTR_mfMCG_iRm=max(mean_env_avg_SUBTR_mfMCG_iRm);

att_mMCG_SUBTR2vsRaw=100-((max_env_avg_SUBTR_fMCG_iRm/max_env_avg_data_iRm)*100);
att_mMCG_SUBTRmfvsRaw=100-((max_env_avg_SUBTR_mfMCG_iRm/max_env_avg_data_iRm)*100);

[yupper,ylower]=envelope(avg_data_iRf);
mean_env_avg_Raw_iRf=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_Raw_iRf(mm)=mean(yupper(mm,:));
end
max_env_avg_Raw_iRf=max(mean_env_avg_Raw_iRf); 

[yupper,ylower]=envelope(avg_SUBTR_fMCG_iRf);
mean_env_avg_SUBTR_fMCG_iRf=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_SUBTR_fMCG_iRf(mm)=mean(yupper(mm,:));
end
max_env_avg_SUBTR_fMCG_iRf=max(mean_env_avg_SUBTR_fMCG_iRf); 

[yupper,ylower]=envelope(avg_SUBTR_mfMCG_iRf);
mean_env_avg_SUBTR_mfMCG_iRf=ones(size(yupper,1),1);
for mm=1:size(yupper,1)
    mean_env_avg_SUBTR_mfMCG_iRf(mm)=mean(yupper(mm,:));
end

max_env_avg_SUBTR_mfMCG_iRf=max(mean_env_avg_SUBTR_mfMCG_iRf);

att_fMCG_SUBTR2vsRaw=100-((max_env_avg_SUBTR_fMCG_iRf/max_env_avg_Raw_iRf)*100);
att_fMCG_SUBTRmfvsRaw=100-((max_env_avg_SUBTR_mfMCG_iRf/max_env_avg_Raw_iRf)*100);

display('attenuations calculated');

% read marker file
[mrk]=getMarkers(fi);
display('truncating stimulus markers');
% if trigger exists, truncate, create itimesTrunc
for ii=1:length(mrk.itimes)
    % truncate the marker times [stTrunc: endTrunc]
    temp_st=find(mrk.itimes{ii}<=stTrunc);
    temp_end=find(mrk.itimes{ii}>=endTrunc);
    temp=mrk.itimes{ii};
    temp([temp_st,temp_end])=[];
    mrk.itimesTrunc{ii}=temp-stTrunc;
    clear temp;
end

% Amplitude sensitive threshold detection 
% Reject trials with amplitude higher than 2pT(2e-12)
display('Rejecting trials with great amplitude');
existtr17=strncmp(mrk.name,'Tr17',4);
%use the marker that has more triggers
if length(existtr17)>=2 && (existtr17(2)==1)
    nTr=2;
else
    nTr=1;
end
tTr=[0.2 0.8];
itimesTruncnewSUBTR=[];
itimesTruncnewSUBTRmf=[];

for ii=1:length(mrk.itimesTrunc{1,nTr})
    Trnow = mrk.itimesTrunc{1,nTr}(ii);
    if Trnow+round(tTr(2)*sf)>size(SUBTR_fMCG,1)
        break;
    end
    if max(max(abs(SUBTR_fMCG(Trnow-round(tTr(1)*sf):Trnow+round(tTr(2)*sf),ExcCh))))<2e-12
        itimesTruncnewSUBTR=[itimesTruncnewSUBTR Trnow];
    end
    if max(max(abs(SUBTR_mfMCG(Trnow-round(tTr(1)*sf):Trnow+round(tTr(2)*sf),:))))<2e-12
        itimesTruncnewSUBTRmf=[itimesTruncnewSUBTRmf Trnow];
    end
end
mrk.itimesTruncnewSUBTR=itimesTruncnewSUBTR;
mrk.itimesTruncnewSUBTRmf=itimesTruncnewSUBTRmf;

if length(itimesTruncnewSUBTR)>0
    [avg_SUBTR_fMCG_tr17,avg_SUBTR_fMCG_tr17_span]=time_locked_avg(SUBTR_fMCG(:,ExcCh),mrk.itimesTruncnewSUBTR,[-round(tTr(1)*sf):round(tTr(2)*sf)]);
    beginpre= 1;
    endpre = find(avg_SUBTR_fMCG_tr17_span==0);
    avg_SUBTR_fMCG_tr17_bnpf=removeDC(band_pass_filter(0.5,10,avg_SUBTR_fMCG_tr17,sf)); 
    avg_SUBTR_fMCG_tr17_bnpf_bc = avg_SUBTR_fMCG_tr17_bnpf';
    [avg_SUBTR_fMCG_tr17_bnpf_bc] = ft_preproc_baselinecorrect(avg_SUBTR_fMCG_tr17_bnpf_bc, beginpre, endpre);
    avg_SUBTR_fMCG_tr17_bnpf_bc = avg_SUBTR_fMCG_tr17_bnpf_bc';
    [yupper,ylower]=envelope(avg_SUBTR_fMCG_tr17_bnpf_bc(beginpre:endpre,:));
    mean_env_avg_SUBTR_n=ones(size(yupper,1),1);
    for mm=1:size(yupper,1)
        mean_env_avg_SUBTR_n(mm)=mean(yupper(mm,:));
    end
    max_env_avg_SUBTR_n=max(mean_env_avg_SUBTR_n);
    [maxval,indx]=max(max(abs((avg_SUBTR_fMCG_tr17_bnpf_bc(avg_SUBTR_fMCG_tr17_span>round(0.2*sf)&avg_SUBTR_fMCG_tr17_span<round(0.5*sf),:))),[],1));
[samp_no , ch_no] = ind2sub(size(avg_SUBTR_fMCG_tr17_bnpf_bc),find(abs(avg_SUBTR_fMCG_tr17_bnpf_bc)==maxval));
latSUBTR=samp_no;
latencySUBTR=latSUBTR/size(avg_SUBTR_fMCG_tr17_bnpf_bc,1)-0.2;
SNRsubtr=maxval/max_env_avg_SUBTR_n;
end

if length(itimesTruncnewSUBTRmf)>0
    [avg_SUBTR_mfMCG_tr17,avg_SUBTR_mfMCG_tr17_span]=time_locked_avg(SUBTR_mfMCG(:,:),mrk.itimesTruncnewSUBTRmf,[-round(tTr(1)*sf):round(tTr(2)*sf)]);
    avg_SUBTR_mfMCG_tr17_bnpf=removeDC(band_pass_filter(0.5,10,avg_SUBTR_mfMCG_tr17,sf)); 
    avg_SUBTR_mfMCG_tr17_bnpf_bc = avg_SUBTR_mfMCG_tr17_bnpf';
    [avg_SUBTR_mfMCG_tr17_bnpf_bc] = ft_preproc_baselinecorrect(avg_SUBTR_mfMCG_tr17_bnpf_bc, beginpre, endpre);
    avg_SUBTR_mfMCG_tr17_bnpf_bc = avg_SUBTR_mfMCG_tr17_bnpf_bc';     
    [yupper,ylower]=envelope(avg_SUBTR_mfMCG_tr17_bnpf_bc(beginpre:endpre,:));
    mean_env_avg_SUBTRmf_n=ones(size(yupper,1),1);
    for mm=1:size(yupper,1)
        mean_env_avg_SUBTRmf_n(mm)=mean(yupper(mm,:));
    end
    max_env_avg_SUBTRmf_n=max(mean_env_avg_SUBTRmf_n);    
[maxval,indx]=max(max(abs((avg_SUBTR_mfMCG_tr17_bnpf_bc(avg_SUBTR_mfMCG_tr17_span>round(0.2*sf)&avg_SUBTR_mfMCG_tr17_span<round(0.5*sf),:))),[],1));
[samp_no , ~] = ind2sub(size(avg_SUBTR_mfMCG_tr17_bnpf_bc),find(abs(avg_SUBTR_mfMCG_tr17_bnpf_bc)==maxval));
latSUBTRmf=samp_no;
latencySUBTRmf=latSUBTRmf/size(avg_SUBTR_mfMCG_tr17_bnpf_bc,1)-0.2;
SNRsubtrmf=maxval/max_env_avg_SUBTRmf_n;
end

display('ER calculated');
