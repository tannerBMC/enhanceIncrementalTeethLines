%
% [normI,minVal,maxVal] = normalizeImagePercentile(I,p,newMax)
%
% normalize image from [p,100-p] percentile to [0,newMax] range

function [normI,minVal,maxVal] = normalizeImagePercentile(I,p,newMax)

minVal = prctile(single(I(:)),p);
maxVal = prctile(single(I(:)),100-p);

disp(['normalized from [' num2str(minVal,'%.2f') ',' num2str(maxVal,'%.2f') '] to [0,' num2str(newMax) ']'])

normI=newMax*(single(I)-minVal)/(maxVal-minVal);
normI(normI(:)<0)=0;
normI(normI(:)>newMax)=newMax;

