function [Signal_Total_cut]=subtractMethod_optimized(Reference_channel,Signal_channel,Fs,fo)
% [Signal_Total_cut]=subtractMethod_optimized2(Reference_channel,Signal_channel,Fs,fo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subtract method:
%Reference_channel = [N sensors x samples] mat file
%Signal_channel=[N sensors x samples] mat file
%Fs= sampling Freq
%fo=cut frequency of the highpass filter (defaul=0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes from original C code:
% 9/24/2014 @Czaragoza and Descalona: changed in 'makeCosWind' function:
% >>> j=segm:1:points-1 for j=points-segm:1:points-1;
% 10/02/2014 @Czaragoza and Descalona: Gains H are not calculated correctly 
% when more than 3 reference channels due to sum of eq. 2.21 (C code use pointers) 
% and dot product can not be used in Matlab or maybe with the "conj" function 
% eq.2.21 : H(i,y)=L(i,y)-sum_(j=i+1)^(q)[L(i,j)*H(j,y)], then change the code as follows:
% >>> gH=zeros(1,length(rgL));
% >>> gH=complex(gH);
% >>> newgH=rgL(indx,indx)-(sum(rgL(indx,:).*gH));
% 10/03/2014 @Czaragoza and Descalona: Removed few "for" loops
% 12/2/2015 descalona: 
% -removing for loops when #input channels=1 #refs=12; for 10min data
% - Matriz multiplication to obtain the Gains, it was avoid the 'while'
% loop of original code (see page 266 of Random Book) 
% 12/3/2015 @descalona removing for loops and split into function calls 
% #input channels=1 #refs=12; for 10min data 
% 12/8/2015 @descalona
% change order of for loops in line >>refSum(g,s)= refSum(g,s)+(Yo(orderingMatrix(s,i),g,s)*gainsH(s,i));
% >>(TW,offs,G)=(780,390,483)
% Elapsed time is 36.020489 seconds.
% >>(TW,offs,G)=(780,390,483)
% Elapsed time is 8.922388 seconds.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shaping = 50;                   % cosine window shaping fixed at 50%
shapeCorr = 3/8;                % fft correction for shaping of 50%
numDeshapZeros = 2;             % number of zeros to add when trail joining
L=size(Signal_channel,2);       % length data

overLap = 0.5;                  % overlap 50%

% frequency resolution (fres) and cut frequency of the high passfilter (fo)
% should be fres<=fo
if nargin>3
   fres=fo;                   
else
   fres = 0.4;   
end
    
TW = round(Fs/fres);             % number of points for fft 
if (mod(TW,2))~=0;TW=TW-1;end   % make even
                
noverlap=ceil(TW*overLap);      % number of sample for overlap
if (mod(noverlap,2)) ~= 0 
    noverlap = noverlap-1;          % make even
end

offset = 1:noverlap:L-TW;      % offset vector 

Nref=size(Reference_channel,1);     % number of reference channels
Nsignal=size(Signal_channel,1);     % number of channels to clean
RC=Nref+1;                          % reference channels + each channels to clean
G=length(offset);                   % number of trials
[window,~]=makeCosWind(shaping,TW); % Create the data segment shaping window

str=sprintf('>>(TW,offs,G)=(%d,%d,%d)',TW,noverlap,G); 
disp(str);

b = zeros (RC,fix(G),TW);
Yo=complex(b);clear b;
for i=1:Nref
    for g=1:G
        %Reference channels trials 
        Ref_trials= Reference_channel(i,offset(g):offset(g)+TW-1);
        %remove the mean
        Ref_trials=Ref_trials-mean(Ref_trials);
        % appply shaping window to each trial
        Out=window.*Ref_trials; 
        % Fourier Transform for each trial and normalize by TW=size of
        % window
        Yotemp=(fft(Out,TW))/TW;
%         plotPSD(Yotemp./sqrt(shapeCorr),Fs,1)    
        % shaping correction
        Yo(i,g,:)=Yotemp./sqrt(shapeCorr);
        clear Out Yotemp
    end
end


for c = 1:1:Nsignal  
    for g=1:G
        % Signal trials
        Signal_trials=Signal_channel(c,offset(g):offset(g)+TW-1);
        %remove the mean
        Signal_trials=Signal_trials-mean(Signal_trials);
        % appply shaping window to each trial
        In=window.*Signal_trials;           
        % Fourier Transform for each signal trial and normalize by TW=size of
        % window
        Yotemp=(fft(In,TW))./TW;
%         plotPSD(Yotemp./sqrt(shapeCorr),Fs,2)
        % shaping correction in signal trials
        Yo(RC,g,:)=Yotemp./sqrt(shapeCorr);
        clear In Yotemp
    end
    
    [XX,orderingMatrix]=Sort_Cohe_CrossSpectral(Yo,G,RC,TW,Nref);  
        
    for s=1:1:TW %L
        mtrx= XX(:,:,s);
        Gxx=mtrx(1:Nref,1:Nref);
        Gxy=mtrx(1:Nref,RC);
        warning off;
        gH=zeros(1,RC);
        gH=complex(gH);
        % here see page 266 of Random Data Book, it was avoid 'while' loop
        % for calculated the Gains, the precision of the calculation of the 
        % gains are slightly modified because the inverse of the matriz Gxx
        gH(:,1:Nref)=(Gxx)\Gxy;
        % gains
        gainsH(s,:)=gH;
    end % done conditioning and GainH

 
    % Reorder ref FTTs using orderingMatrix(NFFT,RC)
    % Mult. ref ffts by gains, sum all ref ffts for each time window (second term of eq 2.2 jiri Notes), and subtr from channel ffts
    b=zeros(fix(G),TW);%L
    refSum=complex(b); % initialize the sum
    
        
    for s=1:TW %L
        for i=1:Nref
            for g=1:fix(G)
                refSum(g,s)= refSum(g,s)+(Yo(orderingMatrix(s,i),g,s)*gainsH(s,i));
            end
        end
    end
 
    clear Yo_temp
    b = zeros (RC,fix(G),TW);
    Yo_temp=complex(b);
    for g=1:fix(G)
        for s=1:TW %L
            Yo_temp(RC,g,s)= (Yo(RC,g,s))-(refSum(g,s)) ; 
        end
        
%         temp=Yo_temp(RC,g,:);temp=temp(:);
%         plotPSD(temp,Fs,10,'r') ;
% 
%         temp2=Yo(RC,g,:);temp2=temp2(:);
%         plotPSD(temp2,Fs,10,'b') ;
    end
           

    % deshaping window 
    deWindow=(1./window).*sqrt(shapeCorr);
    deWindow([1:numDeshapZeros, TW-numDeshapZeros:TW])=0;
    
    % invFFT
    for g=1:fix(G)        
         In(g,:) = Yo_temp(RC,g,:); 
%         temp=arragedFFT(In(g,:));
        invFft(g,:)=(ifft(In(g,:),TW)).*TW.*deWindow; %L,TW,'symmetric'
    end
  
    
    Signal_Total(c,:)=join_series(invFft,G,TW,L,noverlap); 
    
%     Signal_Total(c,:)= Signal;

end  % end c for all primary channels to clean

[stTrunc,endTrunc]=TruncateData(size(Signal_Total,2),Fs,fres);

Signal_Total_cut=zeros(size(Signal_Total));
Signal_Total_cut(:,stTrunc:endTrunc)=Signal_Total(:,stTrunc:endTrunc);

