function [condtensor, maxCond] = sb_assiTensorCond_check_hist(mask,nodes,elem,condcell)
% modified version which excludes outliers from the wm conductivity
% mantonak, February 6, 2018

condtensor = zeros(9,size(elem,1));
count = 0;
countMax =0;
for i = 1 : size(elem)
    pos_o = nodes(elem(i,3),:); 
    pos = round(pos_o);
    if (~(mask.anatomy(pos(1),pos(2),pos(3)) == 0))
        count = count+1;
        S = svd(condcell{pos(1),pos(2),pos(3)});
        maxEig(count) = max(S);
        minEig(count) = min(S);
        prodEig(count) = prod(S);
        meanCond(count) = mean(mean(condcell{pos(1),pos(2),pos(3)}));
        maxCond(count) = max(max(condcell{pos(1),pos(2),pos(3)}));
%         if maxCond(count) < 4 %when plotting the maxCond, I noticed 76 outliers and saw graphically a threshold
            for j = 1 : 3
                for k = 1 : 3
                    condtensor((j-1)*3 + k,i) = condcell{pos(1),pos(2),pos(3)}(j,k);
                end
            end
        else
            countMax = countMax + 1;
%         end
    end
end
count
countMax        
% figure, plot(maxCond)
nbin = 90;
figure,
subplot 231
h = histogram(maxCond);
title('histogram of the maximum of conductivity tensor')
subplot 232
h1 = histogram(maxEig);
title('histogram of the maximum eigen value of conductivity tensor')
subplot 233
h2 = histogram(minEig);
title('histogram of the minimum eigen value of conductivity tensor')
subplot 234
h3 = histogram(prodEig);
subplot 235
h4 = histogram(meanCond);
title('histogram of the mean value of conductivity tensor')
mean(meanCond)
end