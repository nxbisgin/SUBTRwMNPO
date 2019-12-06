function [D,R,sf,LabelSens,pathStr,fileName]=getData(fi)
% fi: folder path of the .ds File
% modified 2012 --  Remove trailing zeros on the end of a data set


    [pathStr,fileName] = fileparts(fi);
    % nice filename
    pos=strfind(pathStr,'\');str=pathStr(pos(end-1)+1:end);str(str=='\')='_';
    disp(strcat('reading data...',str));
    ctf=ctf_read(fi);% The data are loaded from file;
    dall=ctf.data;
  
    e=elec_pos();e=fliplr(e);
    sf=ctf.setup.sample_rate; 
    for i=1:size(dall,2)      
        dd=dall(:,i,:);
        dd=squeeze(dd);
        dd=dd(:);
        d(:,i)=dd;clear dd;
    end   
    f=ctf.sensor.label;f=f';
    m=0;
    for i=1:size(e,1)
        for j=1:size(e,2)
            for k=1:length(f)
                if strncmp(f(k,:),e{i,j},4)~=0;
                    m=m+1;
                    D(:,m)=d(:,k); 
                    R(m,1) = ctf.sensor.info(k).coil(1).position.x;
                    R(m,2) = ctf.sensor.info(k).coil(1).position.y;
                    R(m,3) = ctf.sensor.info(k).coil(1).position.z;  
                    LabelSens{m,:}=e{i,j};
                end
            end
        end
    end
    
    % Remove zeros at the end of meg data
    rng = find(sum(D(:,:)==0,2)<length(R));
    D = D(rng,:);
    
    % Remove channels that show signature of unlocking
    % Added 7/24/2012
    % Added second threshold 8/31/2012
    mark_bad=[];
    maxdiffmeg = max(abs(diff(D)));
    mark_bad1 = maxdiffmeg > (10*median(maxdiffmeg));  % old threshold
    mark_bad2 = maxdiffmeg > (median(maxdiffmeg)+5*std(maxdiffmeg));
    mark_bad = or(mark_bad1, mark_bad2);
    if sum(mark_bad)
        disp('Channels with large transients removed')
        disp(LabelSens(mark_bad,:))
        D(:,mark_bad) = [];
        LabelSens(mark_bad,:) = [];
        R(mark_bad,:) = [];
    end


end