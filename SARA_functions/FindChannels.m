function indxInput=FindChannels(channel_labels,varargin)
% find channels in 'channel_labels' array
% find input channels
if nargin>1  
    
    m=0;
    for ch=1:length(channel_labels);
        for ninput=1:length(varargin{1,:})
            if (strncmp(varargin{1,:}(ninput,:),channel_labels(ch,:),4)~=0)
                m=m+1;
                indxInput(m,:)=ch;  
            end
        end
    end
    indxInput=unique(indxInput);
else
    indxInput=1:length(channel_labels);   
    indxInput=unique(indxInput);
end

    
