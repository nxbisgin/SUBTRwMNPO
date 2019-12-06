function Signal=join_series(invFft,G,TW,L,noverlap) 
i=0;
    for g = 1:fix(G) % reconstruct time series      
        if g==1
            for s=1:(TW-((TW-noverlap)/2)) % join first trial after dropping points from end
                i=i+1;              
                Signal(i)= invFft(g,s);                
            end
        elseif (g>1 && g<fix(G)) % join the intermediate trials after dropping points from each end
            for s=((TW-noverlap)/2)+1:((TW-(TW-noverlap)/2))
                i=i+1;
                Signal(i)=invFft(g,s);
               
            end
        elseif (g==fix(G))   % join last trial after dropping points from start
            for s=(TW-noverlap)/2+1: TW
                i=i+1;
                Signal(i)=invFft(g,s);
            end
        end
        
    end
    
      TR = i;
%      % Pad the end of the data
     for i=TR:L
         Signal(i)=0;
     end
   
end