function [indxbests_mMCG,AmMCG]=SelectSens_mMCG_Based_MaxRel(data,sf,channel_labels,ixExclusionSensors,nFeatures)
% select mMCG sensors based on MUTUAL information
% data: [samples X channels]
% sf: sample frequency
% channel_labels: channels names
% ixExclusionSensors: indexes of sensors around head coil (exclusion area)
% nFeatures: number of sensors want to chose, defauld 12 sensors
% @Descalona 2014
    % defauld nFeatures to chose
    if nargin<5
        nFeatures=12;
    else
        nFeatures=nFeatures;
    end
    addpath ./scripts
    % subsample list where is found the mMCG signal
    % to chose at least the follow sensors if Mutual info is shared
    defauldList_mMCG={'MRQ3';'MCQ0';'MLQ3';'MRP1';'MRP2';'MLP3';  'MCM0'};
    % additional list 
    mMCG_subList={'MLL3';'MRL3';'MLN5';'MLN3';'MRN3';'MRN5'; 'MCO0'};
    
    % index for mMCG
    ixbasedTemplate_mMCG=FindChannels(channel_labels,[defauldList_mMCG;mMCG_subList]);
    ixdefauld_mMCG=FindChannels(channel_labels,defauldList_mMCG);
    ixdefauld_mMCG(ismember(ixdefauld_mMCG,ixExclusionSensors))=[];
    
    % exclude from mMEG template
    ixbasedTemplate_mMCG(ismember(ixbasedTemplate_mMCG,ixExclusionSensors))=[];
    % filter data 0.5-40 Hz
    data_bndf = band_pass_filter(0.5,40,data,sf);
    % 2 minutes of data to calculated Mutual Info
    if size(data,1)>round(60*2*sf)
        block=round(60*2*sf);
    else
        % for short data
        block=size(data,1);
    end
    % Template data mMCG
    TemplateData_mMCG=data_bndf(1:block,ixbasedTemplate_mMCG);
    % selected indexes channels that are not in the template
    ix=ones(size(data_bndf,2),1);
    % avoid features of the tamplate sensors for Mutual information 
    ix(ixbasedTemplate_mMCG)=0; 
    % but take features of the defauld list sensors
    ix(ixdefauld_mMCG)=1;
    % avoid exclusion area to don't chose any sensor
    ix(ixExclusionSensors)=0;
    ix=find(ix);
    % selected features based on mRMR for mMCG
    [selectedFeatures_mMCG,classMI_mMCG,TotalScoreFeatures_mMCG]= MaxRel(nFeatures,data_bndf(1:block,ix),TemplateData_mMCG);
    % indexes of the selected features mMCG
    indxbests_mMCG=ix(selectedFeatures_mMCG);
    % save settings in A struc
    % mMCG
    AmMCG.TemplateData=TemplateData_mMCG;
    AmMCG.ixbasedTemplate=ixbasedTemplate_mMCG;
    AmMCG.selectedFeatures=selectedFeatures_mMCG;
    AmMCG.classMI=classMI_mMCG;
    AmMCG.ixdefauld_mMCG=ixdefauld_mMCG;
    AmMCG.TotalScoreFeatures=TotalScoreFeatures_mMCG;
    AmMCG.ix=ix;
    AmMCG.nFeatures=nFeatures;
        
return

% figure;hold on;
% for jj=1:length(channel_labels)
%     channel_labels2(jj,:)=channel_labels{jj}(2:end);
% end
% plot(R(:,1),R(:,2),'.','Color',[0.7 0.7 0.7]);
% plot(R(ixbasedTemplate_mMCG,1),R(ixbasedTemplate_mMCG,2),'r*');
% plot(R(basedTemplate_fMCG,1),R(basedTemplate_fMCG,2),'b*');
% putSensLabels(R(:,1:2),channel_labels2);
% xlabel('cm');ylabel('cm');
% set(findall(gcf,'type','text'),'FontSize',9,'FontName','times');


