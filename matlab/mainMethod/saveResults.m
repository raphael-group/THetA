%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes the results from running the heterogeneity code and
% saves the results to a text file that can be used later to plot.
%
% INPUT:
% fileName - base name for where reusults will be saved
% bestC - cell array of the best C matrices
% bestMu - cell array of best mu values
% bestP - cell array of corresponding simplex points
% allLike - array of likelihood of all solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = saveResults(fileName, bestC, bestMu, bestP, allLike)

lines = cell(1+size(bestC,2)); % plus one for header line
lines{1} = sprintf('#NLL\tmu\tCounts\tp*');

% Iterate through solutions
for i=1:size(bestC,2)
   curC = bestC{i};
   curMu = bestMu{i};
   curP = bestP{i};
   curLikelihood = allLike(i);
   
   % Convert all pieces to the correct type of string
   
   mu_str = getVectorString(curMu);
   
   curC(:,1) = []; %remove normal column
   c_str = getMatrixString(curC);   
   
   p_str = getVectorString(curP);
   
   
   lines{i+1} = sprintf('%f\t%s\t%s\t%s',curLikelihood,mu_str,c_str,p_str);
   
end

outputFile = strcat(fileName,'.results');
fid = fopen(outputFile,'w');

for i=1:length(lines)
   fprintf(fid,'%s\n',lines{i}); 
end

fclose(fid);

end % end main function


function X_str = getVectorString(curX)
   X_str='';
   for j=1:size(curX,1)-1
       X_str = strcat(X_str,num2str(curX(j,1)),',');
   end
   X_str = strcat(X_str,num2str(curX(end,1)));
end % end of getVectorString


function X_str = getMatrixString(curX)
   X_str='';
   for j=1:size(curX,1)-1
      for k=1:size(curX,2) -1
         X_str = strcat(X_str,num2str(curX(j,k)),','); 
      end
      X_str = strcat(X_str,num2str(curX(j,end)),':');
   end
   
   for k=1:size(curX,2) -1
       X_str = strcat(X_str,num2str(curX(end,k)),',');
   end
   X_str = strcat(X_str,num2str(curX(end,end)));
end
