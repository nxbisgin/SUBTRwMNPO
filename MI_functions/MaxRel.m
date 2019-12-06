function [selectedFeatures,classMI,TotalScoreFeatures] = MaxRel(k, featureMatrix, classColumn)
% function [selectedFeatures,classMI,TotalScoreFeatures] = MaxRel(k, featureMatrix, classColumn)

% @DEscalona 2014: 


noOfFeatures = size(featureMatrix,2);
unselectedFeatures = ones(noOfFeatures,1);

classMI = zeros(noOfFeatures,1);
answerFeatures = zeros(k,1);
highestMI = 0;
highestMICounter = 0;
currentHighestFeature = 0;

featureMIMatrix = -(ones(k,noOfFeatures));

%setup the mi against the class
for n = 1 : noOfFeatures
  classMI(n) = mi(featureMatrix(:,n),classColumn);
	if classMI(n) > highestMI
		highestMI = classMI(n);
		highestMICounter = n;
	end
end

answerFeatures(1) = highestMICounter;
unselectedFeatures(highestMICounter) = 0;
TotalScoreFeatures=zeros(noOfFeatures,k);
%iterate over the number of features to select
for i = 2:k
  score = -100;
	currentHighestFeature = 0;
	iMinus = i-1;
  for j = 1 : noOfFeatures
    if unselectedFeatures(j) == 1
        currentScore = classMI(j);
        TotalScoreFeatures(j,i)=currentScore;
        if (currentScore > score)
            score = currentScore;
            currentHighestFeature = j;
        end
    end
 end
	
  if score < 0 
    disp(['at selection ' int2str(j) ' mRMRD is negative with value ' num2str(score)]);
  end

	%now highest feature is selected in currentHighestFeature
	%store it
	unselectedFeatures(currentHighestFeature) = 0;
	answerFeatures(i) = currentHighestFeature;
end

selectedFeatures = answerFeatures;
