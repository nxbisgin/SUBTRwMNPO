function [indxbests_fMCG,AfMCG]=SelectSens_fMCG_Based_mRMRfeaturesfromExclusionArea(data,sf,channel_labels,ixExclusionSensors,indxbests_mMCG,nFeaturesUser)
% select fMCG sensors based on MUTUAL information using the features form
% exclusion area
% data: [samples X channels]
% sf: sample frequency
% channel_labels: channels names
% ixExclusionSensors: indexes of sensors around head coil (exclusion area)
% nFeatures: number of sensors want to chose, defauld 12 sensors
% @Descalona 2014
% defauld nFeatures to chose
    
    
    if nargin<6
        nfeaturesDefauld=25;
        nFeatures=length(ixExclusionSensors)+nfeaturesDefauld;
        flagnFeatures=0;
    else
        nFeatures=length(ixExclusionSensors)+nFeaturesUser;
        flagnFeatures=1;
    end
    
    % to chose at least the follow subsample list of sensors if Mutual info is shared
    defauldList_fMCG={'MRK4';'MLK4';'MRJ1';'MRI4';'MLH3';'MRG2';'MRF1'};
    % additional subsample list 
    fMCG_subList={
        'MCK0';
        'MRJ3';'MLJ3'
        'MCI0';
        'MRH6';'MRH3';'MLH6';
        'MCG0';
        'MRF3';'MLF3';
        'MCE0';
        'MRD3';'MLD3';
        'MCC0';
        'MRA2';'MLA2'};
    % index for fMCG
    ixbasedTemplate_fMCG=FindChannels(channel_labels,[defauldList_fMCG;fMCG_subList]);
    ixdefauld_fMCG=FindChannels(channel_labels,defauldList_fMCG);
    ixdefauld_fMCG(ismember(ixdefauld_fMCG,indxbests_mMCG))=[];    
    % exclude mMCG ref already chose before
    ixbasedTemplate_fMCG(ismember(ixbasedTemplate_fMCG,indxbests_mMCG))=[];
    % filter data 0.5-40 Hz
    data_bndf = band_pass_filter(1,40,data,sf);
    % 2 minutes of data
    block=round(60*2*sf);
    if block<size(data_bndf,1)
        block=block;
    else
        block=size(data_bndf,1);
    end
    % Template data fMCG include exclusion area
    TemplateData_fMCG=data_bndf(1:block,ixbasedTemplate_fMCG);
    % selected indexes channels that are not in the template ix=146
    ix=ones(size(data_bndf,2),1); 
    % turn off mMCG ref channels which where selected previuos step
    % ix=146-12
    ix(indxbests_mMCG)=0;
    % don't use information from the template ix=146-12-23
    ix(ixbasedTemplate_fMCG)=0;
    % turn on subsample list (defauld features) 146-12-23+7 (repeat ixdefauld_fMCG see line-37 )
    ix(ixdefauld_fMCG)=1;
    % sensors that are turn on  146-12-23+7=118 ix
    ix=find(ix);    
    % selected features based on mRMR for fMCG
    [selectedFeatures_fMCG,classMI_fMCG,TotalScore_fMCG,I_cxj_fMCG,Sum_xi_xj_fMCG,Sum_xi_xjFeatures_fMCG,TotalScoreFeatures_fMCG]= mRMR_D(nFeatures,data_bndf(1:block,ix),TemplateData_fMCG);
    % indexes of the selected features fMCG
    indxbests_fMCG=ix(selectedFeatures_fMCG);
    % removed features from exclusion area
    indxbests_fMCG(ismember(indxbests_fMCG,ixExclusionSensors))=[];
    % select the number of user desired reference which are the most relevant and less redundat
    if flagnFeatures
        indxbests_fMCG=indxbests_fMCG(1:nFeaturesUser);
    end
    
    % save settings in A struc
    % fMCG
    AfMCG.TemplateData=TemplateData_fMCG;
    AfMCG.ixbasedTemplate=ixbasedTemplate_fMCG;
    AfMCG.selectedFeatures=selectedFeatures_fMCG;
    AfMCG.classMI=classMI_fMCG;
    AfMCG.TotalScore=TotalScore_fMCG;
    AfMCG.I_cxj=I_cxj_fMCG;
    AfMCG.Sum_xi_xj=Sum_xi_xj_fMCG;
    AfMCG.Sum_xi_xjFeatures=Sum_xi_xjFeatures_fMCG;
    AfMCG.TotalScoreFeatures=TotalScoreFeatures_fMCG;
    AfMCG.ix=ix;
    AfMCG.nFeatures=nFeatures;
    
return