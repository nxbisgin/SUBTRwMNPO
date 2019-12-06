function [mrk]=getMarkers(fi)
% get triggers markers
% fi: folder path of the .ds File
% mrk: Triggers 
% mrk.itimes: time samples

%%% First get data info 
ctf = readCTFds(fi);
m=0;
if isfield(ctf,'mrk')
    mrks=ctf.mrk;
    for i=1:length(mrks)
        mrkname=ctf.mrk(i).Name;
        % only read 'tr17' and 'tr18'
        if strncmp('Tr17',mrkname,4) || strncmp('Tr18',mrkname,4)
            m=m+1;
            mrk.name{m}=ctf.mrk(i).Name;
            disp(strcat('reading marker...',mrk.name{i}));
            % continuous time
            mrk.itimes{m}=round(((ctf.mrk(i).trial-1).*60.+ctf.mrk(i).time).*ctf.res4.sample_rate);
        end
    end
else
    mrk.name=[];
    mrk.itimes=[];
    disp('No exist any marker file...');
end