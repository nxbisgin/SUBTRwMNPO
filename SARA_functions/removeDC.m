function D=removeDC(D)
% D data where size(D,2)=num Channels
% remove the mean 
for i=1:size(D,2)
    D(:,i)=D(:,i)-mean(D(:,i));
end