function [window,mflg]=makeCosWind(shaping,points)

% modified: @Diana j=segm:1:points-1 for j=points-segm:1:points-1
if shaping > 50
    display('Shaping > 50. From "makeCosWind", Aborting...');
    mflg = -1;
    return;
end

segm = floor(points*shaping/100.);

i=1;
for j=0:1:segm-1
    window(i)=0.5*(1-cos((pi/segm)*j));
    i = i+1;
end

for j=segm:1:points-segm-1 %1:1:(points-2*segm)
    window(i)=1;
    i=i+1;
end

for j=points-segm:1:points-1  %0:1:segm-1
    window(i)=0.5*(1+cos((pi/segm)*(j-points+segm)));
    i=i+1;
end

mflg=0;

return
