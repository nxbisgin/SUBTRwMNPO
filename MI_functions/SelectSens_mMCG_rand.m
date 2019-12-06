
function [indxbests_mMCG_rand]=SelectSens_mMCG_rand(channel_labels,ixExclusionSensors,nFeatures)

if nargin<3
        nFeatures=12;
    else
        nFeatures=nFeatures;
end
    
rng('shuffle');
e=elec_pos;
e2=e(1:6,:);e2=e2(:);e2(strcmp('',e2))=[];
mMCGindx=FindChannels(channel_labels,e2);


mMCGindx(ismember(mMCGindx,ixExclusionSensors))=[];
ix=randperm(length(mMCGindx));
indxbests_mMCG_rand=mMCGindx(ix(1:nFeatures));



