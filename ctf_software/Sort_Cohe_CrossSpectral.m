function [XX,orderingMatrix]=Sort_Cohe_CrossSpectral(Yo,G,RC,TW,Nref)

% Calculation of the average Cross Spectral Density Matrix
b = zeros (RC,RC,TW);
X=complex(b);clear b;
for s=1:1:TW%L     
    for i=1:1:RC
        X(i,:,s)=sum(bsxfun(@times,conj(Yo(i,:,s)),(Yo(:,:,s))),2)./G;
    end
end

b = zeros (RC,RC,TW); %L 
XX=complex(b);clear b;
orderingMatrix=zeros(TW,RC);
for s = 1:1:TW  %L
    for i=1:1:Nref   
        % calculation of the squared coherance input signal and reference channels (abs(Pxy).^2)./(Pxx.*Pyy); % Cxy
        ordCohs(i)=real((abs(X(i,RC,s)).^2)./(X(i,i,s).*X(RC,RC,s)));
        % coherence by frequency 
        %tempCoh(i,s)=ordCohs(i);            
    end
    % order coherence per frequency descend form
    [~,newOrderingSeq]=sort(ordCohs,'descend');
    newOrderingSeq(RC) = RC;
    % save coherence for each frequency
    orderingMatrix(s,:)=newOrderingSeq;
    % order spectral density matrix according coherence
    for i=1:1:RC       
%         for j=1:1:RC
            XX(i,:,s)= X(newOrderingSeq(i),newOrderingSeq,s);
%                 disp(['in >> ',num2str(i),',',num2str(j),' change to  ', num2str(newOrderingSeq(i)),',',num2str(newOrderingSeq(j)), ])
%         end
    end
end

return