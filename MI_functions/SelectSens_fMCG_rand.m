
function [indxbests_fMCG]=SelectSens_fMCG_rand(channel_labels,ixExclusionSensors,indxbests_mMCG,nFeatures)

if nargin<3
        nFeatures=11;
    else
        nFeatures=nFeatures;
end
rng('shuffle');
e=elec_pos;
e3=e(1:end,:);e3=e3(:);e3(strcmp('',e3))=[];
fMCGindx=FindChannels(channel_labels,e3);
fMCGindx(ismember(fMCGindx,ixExclusionSensors))=[];
fMCGindx(ismember(fMCGindx,indxbests_mMCG))=[];
ix=randperm(length(fMCGindx));
indxbests_fMCG=fMCGindx(ix(1:nFeatures));