function [outputdata]=Umit_Heart_Extract_forMATfilesSARA(datamat,ulbpm,llbpm,writeActasDataEditorFile,updateOrCreateMarkerFile,markerName)

% folderName has to be full path
% ulbpm (upper limit beat per minute)
% for fetal heartbeat is 210 and
% llbpm (lower limit) is 100
% upper limit for maternal data is 110 and lower limit is 40
% markerName Marker name for DataEditor

% writeActasDataEditorFile is used to set if the user wants to write the
% results in a folder called act_ht folder

% updateOrCreateMarkerFile if set to 1 , results will be written as marker
% file. If there is a marker file in the folder, program will rename the
% existing marker file to markerfile_backup.mrk file.
% if there is a Marker with the same name specified with markerName,
% program will remove it first and replace it with the new calculated one.

% after marker values are specified this version also applies localmaxima
% to find the local maximum point close to index point found by using the
% hilbert transform.


% program only considers the sensors starting with letter m

% known issues. If data range is too small, this may affect the threshold
% calculation. Maximum and minimum is calculated by using the 1000th
% maximum and 1000th minimum value of the data and multiplying it with the
% threshold multipliers.
% if most of the values of the rr interval are below or above threshold
% this may cause problem. After the first analysis, analyst should adjust the
% thresholds ulbpm and llbpm parameters.
%folderName='C:\SARAGroupWork\UltrasoundStudy\Patient3\umit03\35w2d\umit03_SPONT-STAN_20080523_35w2d.ds';
%folderName='/server/fmeg-ix2/data2/jessica/pat05_700Hz-mMCG.ds';
%folderName='/server/fmeg-ix2/data2/jessica/pat05_soundwithoutsound-mMCG.ds';
%folderName='/server/fmeg-ix2/data2/jessica/pat12_longtrial-mMCG.ds';
%folderName='/server/fmegli03/data2/pam/hr_pam/analysis-3/hr88/32w/hr88_VER-STAN_20060830_32w-11-mMCG.ds';
%folderName='/server/fmegli03/data2/pam/hr_pam/analysis-3/hr88/32w/hr88_AER-STAN_20060830_32w-11-mMCG.ds';

if (exist('markerName')==0)
    markerName='fMCG2';
end

if (exist('writeActasDataEditorFile')==0)
    writeActasDataEditorFile=0;
end

if (exist('updateOrCreateMarkerFile')==0)
    updateOrCreateMarkerFile=0;
end


thresholdmultiplier=[0.1:0.05:0.9];
% threshold is set based on the maximum and minimum of the data and
% threshold multiplier.

if (exist('ulbpm')==0)
    ulbpm=210; %upper limit
end

if (exist('llbpm')==0)
    llbpm=100; %for number of beats per minute
end
 %%%%%%%%%% datamat struc
folderName=datamat.Comment;
d=datamat.data;
frequency=datamat.sf;
numberofsamples=size(d,1);
F=datamat.good_channel_labels';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pathStr,fileName] = fileparts(folderName);
if (length(pathStr)<3)
    display('Program stopped because of error.');
    display('Please enter the full path information like c://datafolder//datafile.ds');
    return;
end

niceName=fileName;
for i=1:length(niceName)
    if (strncmp(char(95),niceName(i),1))
        niceName(i)=' ';
    end
end



%% load data info, find the intrested channels

%  eval(sprintf('!cp -a %s /opt_prg/SOFTWARE/P_detection/datfolder/hariact/.',tline))


% issave=datamat.issave;


li=length(F);
for i=1:li
    if (strncmp(F(i,:),'M',1)==1)
        minchannel=i;
        break;
    end
end;

for i=li:-1:1
    if (strncmp(F(i,:),'M',1)==1)
        maxchannel=i;
        break;
    end
end;

%% Calculate Mean and power spectrum for Each Channel
% data2=getTrial2(ds,-1,'t');
% 
% firsttime=1;
% for i=1:size(data2,2)
%     dd=data2(:,i,:);dd=squeeze(dd);dd=dd(:);
%     if firsttime==1
%         d=zeros(length(dd),size(data2,2));
%         firsttime=0;
%     end
%     d(:,i)=dd;clear dd;
% end


% meanarray=zeros(numberoftrials,maxchannel-minchannel+1);
%
% for k=1:numberoftrials
%     meanarray(k,:)=mean(d(:,minchannel:maxchannel));
% end


%remove the mean

for j=minchannel:maxchannel
    d(:,j)=d(:,j)-mean(d(:,j));
end


% Filtering
[cf1 cf2]=butter(4,[1/(frequency/2) 60/(frequency/2)]);
%FILTERING OPERATION
mergedarray=filtfilt(cf1,cf2,d(:,minchannel:maxchannel));


firsttime=1;
bestfreqsum=-999999;
best1=0;
bests=[];
freqsums=zeros(1,size(mergedarray,2));

for i=1:size(mergedarray,2)


    [f1,tsum1]=powerspec(mergedarray(:,i), frequency);

    if firsttime==1
        psize=round((90/(frequency/2))*length(tsum1));
        pos60=round((90/(frequency/2))*length(tsum1));
    end

    freqsum=sum(tsum1([1:pos60-10 pos60+10:psize]));

    freqsums(i)=freqsum;

    if freqsum>bestfreqsum
        best1=i;
        bestfreqsum=freqsum;
    end

end

[ffre ind]=sort(freqsums,'descend');

%take the top 10 highest frequency components
bests=ind(1:10);

amplitudearray=d(:,best1+minchannel-1)';

%% after memory allocation


% stores hilbert transformed data

% data is loaded as seperate trials, each load ; 2 trials are loaded and hilbert
% transfor is applied to those trials. Then they are merged. Objective is to get rid of the calculation
% errors of the hilbert transform at the beginning and end of the 2 trial
% series.

% hilbert transformed data is 4th order butterwort filtered
% remove the mean of each trial
%




%FILTERING OPERATION
% Hilbert for polarization
%hilbert transform part
hilbertdata=hilbert(mergedarray(:,bests));
hilbertaxx=abs(diff(hilbertdata));
hilberta=sum(hilbertaxx')';


%% detect max points

%I use the value of 1000 th point and from the end 1000th point
%to determine the min and max average and threshold point


%%%%Length may change


%because of the hilbert transform it may have high values in the beginning
tempaa=sort(hilberta(180:end),'descend');

upperandlowerthreshold=round(length(tempaa)*0.01);
minaverage=tempaa(end-upperandlowerthreshold);
maxaverage=tempaa(upperandlowerthreshold);
clear tempaa;


results=cell(length(thresholdmultiplier),5);
%max points
%amplitude of max points
%thresholod multiplier
%threshold value
%folderName
for thresholdindex=1:length(thresholdmultiplier)

    threshold=minaverage+(maxaverage-minaverage)*thresholdmultiplier(thresholdindex);
    [results{thresholdindex,1}, results{thresholdindex,2}]=getMaximumPointsAboveThreshold(hilberta,threshold);

end

results{thresholdindex,3}=(thresholdmultiplier(thresholdindex));
results{thresholdindex,4}=(threshold);
results{thresholdindex,5}=(folderName);



%% Clear the wrong R peaks
% try to detect the wrong max points and remove them from the array
% if time difference between points are below certain level which causes
% abnormal number of beat per minute.
% missed beat can be because of two things, first there is a small peak
% before the real peak or there can be a wrong counted peak after the real
% peak. Eather way what we can do is to calculate two bear per minute value
% and pick the one closest to the previous number of beats. Because we
% don't expect sudden changes and the next heartbeat value is more likely
% close to the previous heartbeat value.


considerAmplitude=0;


bestthresholdindex=-1;
bestthresholdscore=999999;

for thresholdindex=1:length(thresholdmultiplier)

    dp=diff(results{thresholdindex,1})./frequency;

    dp=60./dp;

    mis1 = sum(dp>ulbpm);
    ext1 = sum(dp<llbpm);
    score=mis1+ext1;

    if (score<bestthresholdscore)
        bestthresholdscore=score;
        bestthresholdindex=thresholdindex;
    end
end

temptemp=1/frequency;
tempend=(length(hilberta)-1)*temptemp;


dp=diff(results{bestthresholdindex,1})./frequency;
meandp=mean(dp);

meanampl=mean(hilberta(results{bestthresholdindex,1})); %mean amplitude




thr1=60./llbpm;
thr2=60./ulbpm;


%if there is anormaly calculate current+next    and next+nextnext to see
%which two values are close to the local average at that region

activeindex=results{bestthresholdindex,1};

if bestthresholdscore==0
    indexnew=activeindex;
else

    for kkk=1:2


        indexnew=[activeindex(1:2)];
        lengthactiveindex=length(activeindex);
        dpnew=dp(1);

        for rty=3:lengthactiveindex


            indexnew=[indexnew activeindex(rty)];
            dpnew=[dpnew dp(rty-1)];

            previous=(indexnew(end-1)-indexnew(end-2))./frequency;
            current=(indexnew(end)-indexnew(end-1))./frequency;;
            indexofprevious=indexnew(end-1);
            indexofcurrent=indexnew(end);

            if (lengthactiveindex==rty)
                indexofnext=length(hilberta);
            else
                indexofnext=activeindex(rty+1);
            end



            if (current>thr1)  && (kkk==1)
                %lets look at the original data and find the maximum point
                %in this current interval

                if (length(indexnew)<11)
                    tempthreshold=meanampl*0.5;
                else
                    tempthreshold=mean(hilberta(indexnew(end-10:end-1)))*0.5; %mean amplitude of the last 10 beats
                end


                rrr=getMaximumPointsAboveThreshold(hilberta(indexofprevious:indexofcurrent),tempthreshold);



                rrr=rrr+(indexofprevious-1);

                timediffs=diff([indexnew(end-1) rrr(2:end) indexofcurrent])./frequency;
                dpnew=[dpnew(1:end-1) timediffs];


                indexnew=[indexnew(1:end-1) rrr(2:end) indexofcurrent];


            end


            if (((current<thr2)||(previous<thr2)) && (kkk>1))

                if (length(dpnew)<12)
                    localmean=meandp;
                else
                    localmean=mean(dpnew(end-11:end-2));
                end


                %we have a small time difference then question is should we add
                %this short time to the previous value or merge it with the
                %next value. One answer might be if we merge it with the
                %previous value then is the value of the previous value becomes
                %closer to the local mean at that region or not
                if ((abs(localmean-(previous+current))<abs((localmean-current)))||(abs(localmean-(previous+current))<abs((localmean-previous))))
                    %then lets merge it because new value of the previous
                    %becomes more close to the local average at that region
                    val1=previous+current;
                    dpnew=dpnew(1:end-1);
                    %remove the last element
                    %update the previous one
                    dpnew(end)=val1;

                    %indexnew(end-1)=indexnew(end); %use the previous index
                    %and update the previous index's value
                    indexnew(end-1)=indexnew(end);
                    indexnew=indexnew(1:end-1); %just remove that index



                    %update previous
                end
            end



        end  % for

        dp=dpnew;
        clear dpnew;
        check=max(activeindex);
        activeindex=indexnew;
        check=max(activeindex);
    end


end %if

indexold=results{bestthresholdindex,1};

time_vector=indexnew./frequency;
  
%% Load highest amplitude channel data and plot yy graph


filteredAmplitudeArray=amplitudearray;%filtfilt(cf1,cf2,amplitudearray);%
sensorindexnew=updateMarkerPointsByUsingLocalMaxima(filteredAmplitudeArray,indexnew);
%figure;hist(sensorindexnew-indexnew,-9:1:9);


diffindexold=diff(indexold)./frequency;
diffindexnew=diff(indexnew)./frequency;


% figure
% subplot(4,1,1)
% plot(indexold(1:end-1)./frequency,(60./diffindexold),'.');
% title(strcat('Original Without Correction - ',niceName))
% %plot(deletedIndexList,hilberta(deletedIndexList),'^');
% ylabel('HeartRate (bpm)');
% xlabel('Time (Second)')
% 
% subplot(4,1,2)
% plot(indexnew(1:end-1)./frequency,(60./diffindexnew),'+');
% title('After Correction with Adaptive Threshold')
% ylabel('HeartRate (bpm)');
% xlabel('Time (Second)')
% 
% 
% subplot(4,1,3)
% plot(indexnew(1:end-1)./frequency,hilberta(indexnew(1:end-1)));
% title('Hilbert Amplitude')
% ylabel('Hilbert Amplitude');
% xlabel('Time (Second)')
% 
% subplot(4,1,4)
% %shift one point back because hilbert and diff is causing one point shift
% plot(sensorindexnew(1:end-1)./frequency,filteredAmplitudeArray(1,(sensorindexnew(1:end-1))),'r');
% ylabel('Sensor Amplitude');
% xlabel('Time (Second)')
% title(strcat('Sensor Measurement -  (',F(best1,:),')'));

% close all

%% Write Results to the Marker File
%calculateTrialandTimes

%keyboard;

if (updateOrCreateMarkerFile==1)

    trials=floor(indexnew./numberofsamples);
    % shift one back because the hilbert transform and diff causes 1 point
    % caused problem with time senc.
    timesec=((indexnew)-(trials.*numberofsamples))./frequency;

    eval(sprintf('!cp %s/MarkerFile.mrk %s/MarkerFile_BACKUP.mrk',folderName,folderName))
    deleteMarkerFromMarkerfile(folderName,markerName);
    writeMarkerfile(folderName,markerName,[trials' timesec'])

    %eval(sprintf('!mkdir /opt_prg/SOFTWARE/P_detection/MarkerFiles/%s',fileName))
    %ctf_write_markerfile(strcat('/opt_prg/SOFTWARE/P_detection/MarkerFiles/',fileName),ctf,'fMCG',[trials' timesec']);
end

%% Write Results
hrvalues=60./diffindexnew;


if (writeActasDataEditorFile==1)

    aa=readDs2(folderName);
    aa.res4.no_channels=3;
    aa.res4.chanNames=[aa.res4.chanNames(51,:);aa.res4.chanNames(51,:);aa.res4.chanNames(51,:)];
    fileNamelength=length(aa.res4.chanNames(1,:))/2;
    for eqw=1:fileNamelength
        aa.res4.chanNames(1,eqw)=' ';
        aa.res4.chanNames(2,eqw)=' ';
    end
    aa.res4.chanNames(1,1)='f'; aa.res4.chanNames(1,2)='H';
    aa.res4.chanNames(1,3)='R';

    aa.res4.chanNames(2,1)='h';
    aa.res4.chanNames(2,2)='a'; aa.res4.chanNames(2,3)='c';

    aa.res4.chanNames(3,1)='s';
    aa.res4.chanNames(3,2)='a'; aa.res4.chanNames(3,3)='c';


    aa.res4.senres=[aa.res4.senres(51) ;aa.res4.senres(51);aa.res4.senres(51)];
    %%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%



    aa.sensor.info(:).proper_gain=1;
    aa.sensor.info(:).q_gain=1;
    aa.res4.senres(1,1).properGain=1;
    aa.res4.senres(1,1).qGain=1;
    aa.res4.senres(1,1).gain=1;

    aa.res4.senres(2,1).properGain=1;
    aa.res4.senres(2,1).qGain=1;
    aa.res4.senres(2,1).gain=1;


    aa.res4.senres(3,1).properGain=1;
    aa.res4.senres(3,1).qGain=1;
    aa.res4.senres(3,1).gain=1;



    %%%%%%%%%%%AMPLITUDE VALUES ARE MULTIPLIED WITH 200./meanamplitudearray
    %%%%%%%%%%%%%%%%%% AND 200./meanhilberta
    %%%%%%%%%%%TO BE STORED AS INTEGERS BECAUSE%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%WRITEDS2 USES INT32 TO STORE THE VALUES%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    meanhilberta=mean(hilberta(indexnew(1:end-1)));
    meanamplitudearray=mean(abs(amplitudearray(1,sensorindexnew(1:end-1))));

    %set the average around 200 for display reasons

    col3=abs((amplitudearray(1,sensorindexnew(1:end-1)).*(200./meanamplitudearray)))';
    %amplitude array


    col2=hilberta(indexnew(1:end-1)).*(200./meanhilberta);
    fileNameofthefhractfile = [fileName '_fhr_act_ht'];
    [cv,vb]=writeDs2(aa.res4,[hrvalues' col2 col3],['int'],fileNameofthefhractfile);


end

%% return values

outputdata.hilberta=hilberta;
% outputdata{2}=[hrvalues' hilberta(indexnew(1:end-1)) (filteredAmplitudeArray(1,sensorindexnew(1:end-1)))'  indexnew(1:end-1)'];
outputdata.indexold=indexold;
outputdata.indexnew=indexnew;
outputdata.amplitudearray=amplitudearray;
outputdata.chanNames=F(best1,:);
outputdata.sensorindexnew=sensorindexnew;
outputdata.folderName=folderName;
outputdata.time_vector=time_vector;

return


function [maxpoints amplitudes]=getMaximumPointsAboveThreshold(datapoints,threshold)
%threshold=minaverage+(maxaverage-minaverage)*thresholdmultiplier(thresholdindex);
maxpoint=1;
maxpoints=[];
localmax=datapoints(1);
stored=localmax<threshold;
for i=1:length(datapoints)
    if (datapoints(i)>threshold)


        if (datapoints(i)>localmax)
            localmax=datapoints(i);
            maxpoint=i;
            stored=0;
        end

    else
        if (stored==0)
            %store the identified max point
            maxpoints=[maxpoints maxpoint];
            stored=1;
        end
        localmax=-999999;
    end %if
end %for

%try to find the point where the values goes below the threshold
%save the local max point and keep progressing until the value goes
%above the treshold

amplitudes=datapoints(maxpoints);
return

function x=gettrialvalues(inarray,i);
%18750 samples for each trial
mynus=(i-1)*18750;
mynusend=i*18750+1;
mystart=find(mynus<inarray ,1,'first');
myend=find(inarray<mynusend,1,'last');

x=inarray(mystart:myend)-mynus;


return

function [f,powe]=powerspec(y, Fs)
%Fs=frequency;
%y=t;
%for sara data y=y'
% y is the data L is window length for fft and Fs samplping frequency
L=length(y);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%,NFFT
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2);
powe=2*abs(Y(1:NFFT/2,:));
%powe=Y(1:NFFT/2);
% Plot single-sided amplitude spectrum.
%length(f)./max(f)*90
return

function [foundMarkerPoints]=updateMarkerPointsByUsingLocalMaxima(sensorData,foundMarkerPoints)
localmaximawindowlength=20;
sensorDataLength=length(sensorData);
for i=1:length(foundMarkerPoints)
    point=foundMarkerPoints(i);

    if ((sensorDataLength>point+10)&&(point>10))
        windowData=sensorData(point-localmaximawindowlength+10:point+10);
        res=find(windowData==max(windowData),1,'last');
        foundMarkerPoints(i)=res+point-localmaximawindowlength+10-1;
    end
end
return

return

function []=deleteMarkerFromMarkerfile(folder,marker_name)
% trials start from 0
% ctf_write_markerfile - write data to MarkerFile.mrk in .ds folder
%
% [] = ctf_write_markerfile([folder],[ctf]);
%
% folder is a .ds path, if it is omitted a gui prompts for the folder.
% ctf is a struct generated by this and other ctf_read_*.m functions.
%
% Returns nothing:
%
% ctf.markers.number_markers - scalar
% ctf.markers.number_samples - column array, number_markers x 1
% ctf.markers.marker_names - cell array of strings, number_markers x 1
% ctf.markers.trial_times - cell array of matrices, number_markers x 1 each
% containing number_samples x 2 matrix, where column 1 indicates the trial
% number containing the marker and column 2 indicates the offset (in sec).
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%

% $Version: 1.00 $ $Date: 2006/03/19 02:39:57 $

% Copyright (C) 2008  Umit D. Ulusar
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
%
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 12/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - modified from NIH code readmarkerfile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Version: 1.00 $';
fprintf('\nCTF_DELETE_MARKER_FROM_MARKERFILE [v %s]\n',ver(11:15)); tic;

% check the inputs
ctf = ctf_folder(folder);
fprintf('...writing to %s\n',ctf.folder);

% check for a marker file
% if we are going to use a seperate marker file function has to take it as
% an input

nomarkerfile=0;
markerFile = findmrkfile( ctf.folder, filesep);
a = exist(markerFile);
if a == 0,
    fprintf('...MarkerFile.mrk does not exist in this .ds folder. ');
    nomarkerfile=1;
    return
end

% read the marker file, delimited by carriage returns.  This new array is
% used below to extract relevant information

if (~nomarkerfile)
    fprintf('...reading all marker file text, using ''\\n'' delimiter\n');
    markerFileText = textread(markerFile,'%s','delimiter','\n');

    % extract the number of markers
    number_markers_index = strmatch('NUMBER OF MARKERS:',markerFileText,'exact');
    number_markers = str2num(markerFileText{number_markers_index + 1});
    fprintf('...found %d marker types:\n',number_markers);

    % extract marker names (this can be an array, if number_markers > 1)
    marker_names_index = strmatch('NAME:',markerFileText,'exact');
    marker_names = markerFileText(marker_names_index + 1);
    for m = 1:length(marker_names),
        fprintf('   %40s\n',marker_names{m});
    end

end %if
%% delete records





newMarkerFileText=deleteMarkerwithMarkerName(markerFileText,marker_name);

if length(newMarkerFileText)>=length(markerFileText)
    % no marker with that name can be found
    return;
end

number_markers = str2num(newMarkerFileText{number_markers_index + 1})
number_markers=number_markers-1
number_markers_index = strmatch('NUMBER OF MARKERS:',newMarkerFileText,'exact');
newMarkerFileText{number_markers_index + 1}=num2str(number_markers);


class_id_index = strmatch('CLASSID:',newMarkerFileText,'exact');
for i=1:length(class_id_index)
    newMarkerFileText{class_id_index(i)+1}=num2str(i);
end





fprintf('Marker is deleted. Now we have %d marker types:\n',number_markers);

%write everything
if (nomarkerfile)
    markerFile=[folder '\MarkerFile.mrk'];
end

fid=fopen(markerFile,'w');


for i=1:length(newMarkerFileText)
    fprintf(fid,'%s\n',newMarkerFileText{i});
end
fclose(fid);




t = toc; fprintf('...done (%6.2f sec)\n\n',t);



return

% -------------------------------------------------------
% function mrkname = findmrkfile( folder, filesep)
% mrkname = dir([ folder filesep 'MarkerFile.mrk' ]);
%
% if isempty(mrkname)
%     fprintf('...no file with extension .mrk or .MRK in .ds folder\n');
%     mrkname = [];
% else
%     mrkname = [ folder filesep mrkname.name ];
% end;
% return

function [newMarkerFileText]=deleteMarkerwithMarkerName(markerFileText,markerName)
marker_names_index = strmatch('NAME:',markerFileText,'exact');
marker_names = markerFileText(marker_names_index + 1);
deleteIndex=-1;
for m = 1:length(marker_names),
    if (strcmp(marker_names{m},markerName)==1)
        deleteIndex=m;
    end
end

if (deleteIndex==-1)
    fprintf('...no Marker with the name %s\n',markerName);
    newMarkerFileText=markerFileText;
    return
end


newMarkerFileText=cell(1,1);
for i=1:(marker_names_index(deleteIndex)-3)
    newMarkerFileText{i}=markerFileText{i};
end

if length(marker_names_index)>=(deleteIndex+1)
    for i=marker_names_index(deleteIndex+1)-2:length(markerFileText)
        newMarkerFileText{end+1}=markerFileText{i};
    end
end
%     newMarkerFileText{1:marker_names_index(deleteIndex)-3}=markerFileText{1:marker_names_index(deleteIndex)-3};


%markerFileText{marker_names_index(deleteIndex)-2:marker_names_index(deleteIndex+1)-3}
return

function [ctf] = ctf_folder(folder,ctf);

% ctf_folder - get and check CTF .ds folder name
%
% [ctf] = ctf_folder( [folder], [ctf] );
%
% folder:  The .ds directory of the dataset.  It should be a complete path
% or given relative to the current working directory (given by pwd).  The
% returned value will ensure the complete path is identified.  If this
% argument is empty or not given, a graphical prompt for the folder
% appears.
%
% eg,
%     ctf = ctf_folder;
%
% ctf.folder is returned (as a complete path).
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%


% $Revision: 1.8 $ $Date: 2005/12/06 05:28:45 $

% Copyright (C) 2004  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('folder','var'), folder = []; end
if ~exist('ctf','var'), ctf = []; end

if isempty(folder),
    if isfield(ctf,'folder'),
        folder = ctf.folder;
    else
        folder = [];
    end
end

if exist(folder) ~= 7,
    fprintf('...folder inputs are invalid\n');
    folder = getfolder;
end

ctf.folder = folder;

% ensure we get the folder path
current_dir = pwd;
cd(ctf.folder);
cd ..
folderPath = pwd;
cd(current_dir);

% check whether the folder path is in the folder already
[path,file] = fileparts(ctf.folder);

% if findstr(folderPath,ctf.folder),
% OK the path is already in the folder
% else
if isempty(path)
    % Add the folderPath
    ctf.folder = fullfile(folderPath,ctf.folder);
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function folder = getfolder,

% This does not work on linux systems (only windows)
folder = uigetdir(pwd,'locate CTF .ds folder');

if ~folder,
    % try to get it this way
    [filename, folder, filterindex] = uigetfile('*.res4', 'Pick a .res4 file');
    if ~folder,
        error('invalid folder');
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function [] = writeMarkerfile(folder,newmarker_name,trial_times)
% trials start from 0
% ctf_write_markerfile - write data to MarkerFile.mrk in .ds folder
%
% [] = ctf_write_markerfile([folder],[ctf]);
%
% folder is a .ds path, if it is omitted a gui prompts for the folder.
% ctf is a struct generated by this and other ctf_read_*.m functions.
%
% Returns nothing:
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%

% $Version: 1.00 $ $Date: 2006/03/19 02:39:57 $

% Copyright (C) 2008  Umit D. Ulusar
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
%
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Version: 1.00 $';
fprintf('\nwriteMarkerFile [v %s]\n',ver(11:15)); tic;

fprintf('...writing to %s\n',folder);

% check for a marker file
% if we are going to use a seperate marker file function has to take it as
% an input

nomarkerfile=0;
markerFile = findmrkfile( folder, filesep);
a = exist(markerFile);
if a == 0,
    fprintf('...MarkerFile.mrk does not exist in this .ds folder. \n New will be created');
    %ctf.markers = [];
    %return
    nomarkerfile=1;
end

% read the marker file, delimited by carriage returns.  This new array is
% used below to extract relevant information

if (~nomarkerfile)
    fprintf('...reading all marker file text, using ''\\n'' delimiter\n');
    markerFileText = textread(markerFile,'%s','delimiter','\n');

    % extract the number of markers
    number_markers_index = strmatch('NUMBER OF MARKERS:',markerFileText,'exact');
    number_markers = str2num(markerFileText{number_markers_index + 1});
    fprintf('...found %d marker types:\n',number_markers);

    % extract marker names (this can be an array, if number_markers > 1)
    marker_names_index = strmatch('NAME:',markerFileText,'exact');
    marker_names = markerFileText(marker_names_index + 1);
    for m = 1:length(marker_names),
        fprintf('   %40s\n',marker_names{m});
    end

end %if
%% add sent records

if nomarkerfile
    initialheader={'PATH OF DATASET:';folder;'';'';'NUMBER OF MARKERS:';'1';'';'';}; %if no marker file than we can use this to initialize
    number_markers=0;
end

extra={'CLASSGROUPID:';'3';'NAME:';char(newmarker_name);'COMMENT:';'Markers of R points';'COLOR:';'red';'EDITABLE:';'Yes';'CLASSID:';num2str(number_markers+1);'NUMBER OF SAMPLES:';num2str(length(trial_times));'LIST OF SAMPLES:';'TRIAL NUMBER		TIME FROM SYNC POINT (in seconds)'}

for i=1:size(trial_times,1)
    extra{end+1}=sprintf('+%s				            +%s',num2str(trial_times(i,1)),num2str(trial_times(i,2)));
end

extra{end+1}=''; %add two blank lines
extra{end+1}=''; %add two blank lines

number_markers=number_markers+1;


if (~nomarkerfile)

    number_markers_index = strmatch('NUMBER OF MARKERS:',markerFileText,'exact');
    markerFileText{number_markers_index + 1}=num2str(number_markers);
    number_markers = str2num(markerFileText{number_markers_index + 1});
end

fprintf('new marker is added. Now we have %d marker types:\n',number_markers);

%write everything
if (nomarkerfile)
    markerFile=[folder filesep 'MarkerFile.mrk'];
end

fid=fopen(markerFile,'w');

if (~nomarkerfile)
    for i=1:length(markerFileText)
        fprintf(fid,'%s\n',markerFileText{i});
    end
else
    for i=1:length(initialheader)
        fprintf(fid,'%s\n',initialheader{i});
    end

end
for i=1:length(extra)
    fprintf(fid,'%s\n',extra{i});
end
fclose(fid);
t = toc; fprintf('...done (%6.2f sec)\n\n',t);


return

% -------------------------------------------------------
% find file name if truncated or with uppercase extension
% added by Arnaud Delorme, June 15, 2004
% edited by Umit Ulusar, April 1 2008 to find only MarkerFile.mrk
function mrkname = findmrkfile( folder, filesep)
mrkname = dir([ folder filesep 'MarkerFile.mrk' ]);
if isempty(mrkname)
    mrkname = dir([ folder filesep 'MarkerFile.MRK' ]);
end;
if isempty(mrkname)
    fprintf('...no file with the name of MarkerFile.mrk or MarkerFile.MRK in .ds folder\n');
    mrkname = [];
else
    mrkname = [ folder filesep mrkname.name ];
end;
return
