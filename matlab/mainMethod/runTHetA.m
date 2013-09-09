%
% 2013 Brown University, Providence, RI.
%
%                       All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of Brown University not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
%
% BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
% http://cs.brown.edu/people/braphael/software.html
% 
% @author Layla Oesper, Ahmad Mahmoody and Benjamin J. Raphael
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes in the following input:
% 
% INPUT:
% inputfile: 5 columns of data: ID, chrm, start, end, tumorCount, normalCount
% inputfile: 7 columns of data (same as above) plus upper and lower bounds
%               on interval counts.
% N: the number of genomes in the mixture to solve for
% K: the max value of k to solve for
% TAU: the expected number of copies in a normal genome
%
% OUTPUT:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[bestLike, bestMu, bestC, bestP, allLike] = runTHetA(inputFile, N, K, TAU, genome, exitVal)

% Step 0: Turn off warnings
s = warning('off','all'); %Turn off warning for performance improvement

% Step 1: add access to all tools necessary (minConf etc.)
file_path = which('runTHetA.m');
f = filesep;
sep_locs = regexp(file_path,f);
idx = sep_locs(end - 1); %Go up one level
tools = strcat(file_path(1:idx),'tools');
addpath(genpath(tools));

% Step 2: make sure input is in correct format
[N,K,TAU,exitVal] = parseInput(N,K,TAU,exitVal);

% Step 3: Load in all data
data = importdata(inputFile);
%If header line was included, we need to get the correct field
if (isfield(data,'data'))
    data = data.data;
end

%Check to see if ID column is still included

%if 6 columns then, remove the first one
if (size(data,2) == 6)
    data = data(:,2:end);
end

%if 8 columns then has bounds and ID, remove first
if (size(data,2) == 8)
    data = data(:,2:end);
end

% Get the variables that we need
M = size(data,1);
tumor = data(:,4);
normal = data(:,5);

num = size(data,2);
if (num>=6)
    bounds = data(:,6);
else
    bounds = K * ones(M,1);
end

if (num==7)
    lowerBounds = data(:,7);
else
    lowerBounds = zeros(M,1);
end

firstCol = TAU * ones(M,1);
init_mu = (1/N) * ones(N,1);

% Values to store when running optimization
bestLike = Inf;
bestMu = {};
bestC = {};
bestP = {};
allLike = [];


epsilon = 0.1;

%initialize first matrix where all columns have 1 in first row
%C_tumor = zeros(M,N-1);
%C_tumor(1,:) = ones(1, N-1);
C_tumor = repmat(lowerBounds,1,N-1);

hasNext = 1;
idx = 1;

%Get sorted information
ratios = tumor./normal;
[a,sorted_idx] = sort(ratios);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Iterate through all valid C and run optimization

tic;

while(hasNext)

    ret = mod(idx,1000);
    if (ret==0)
        fprintf('Trial: %d\n',idx);
    end
     
    C = [firstCol,C_tumor];

    % minFunc options
    opts.verbose = 0;
    opts.numDiff = 0; % use gradient to calculate
    opts.optTol = 1e-4;
    opts.maxIter = 50;
    opts.feasibleInit=1;
    %opts.progTol = 1e-4;
    %opts.testOpt = 0;
 

    
    funObj = @(mu)multinomialLikelihood(mu,C,tumor,normal,N);
    funProj = @(mu)projectSimplex(mu); %restrict solutions to the simplex
    
    %tic;

    [final_mu, likelihood] = minConf_SPG(funObj, init_mu, funProj, opts);
    
    %opts.corrections = 10;
    
    
    %init_mu_lambda = [init_mu;10000];
    %lb = zeros(N+1,1);
    %up = inf(N+1,1);
    
    %funObj2 = @(mu_lambda)laGrangeLikelihood(mu_lambda,C,tumor,normal,N);
    
    %tic;
    %[final_mu_lambda,likelihood] = minConf_TMP(funObj2,init_mu_lambda,lb,up,opts);
    %[final_mu] = doFastOptimization(C,tumor,normal,0,1,0.000001);
    %likelihood = multinomialLikelihood(final_mu,C,tumor,normal,N);  
    
    
    %t1 = toc;
    %fprintf('Time: %s\n',t1);
    
    %final_mu = final_mu_lambda(1:end-1);
    %final_mu = final_mu ./sum(final_mu);
    
    
    
    %t2 = toc;
    %fprintf('Time: %s\n',t2);

    %likelihood = multinomialLikelihood(final_mu,C,tumor,normal,N);
    
    P = getSimplexPt(C,final_mu,normal);
   

    % Store the best solutions
    if (likelihood < bestLike)
       bestLike = likelihood;
    end
    

    if (likelihood <= bestLike + epsilon)

        bestMu{size(bestMu,2) + 1} =final_mu;
        bestC{size(bestC,2) + 1} = C;
        bestP{size(bestP,2) + 1} = P;
        allLike = [allLike; likelihood];
    end

    
    isValid=0;
    while(~isValid && hasNext)
        %Get the next matrix
        [hasNext, C_tumor] = generateNextC(C_tumor,K, bounds, lowerBounds); 
        idx = idx + 1;
        isValid = checkValidC(C_tumor,sorted_idx);
        %isValid=1;
    end
end

%t2 = toc;
%fprintf('Time: %s\n',t2);


%Prune results
[bestC,bestMu,bestP,allLike] = filterBests(bestC,bestMu,bestP,allLike,bestLike, epsilon);


% Step 5: Plot and save all results files

%Plot the Results
fileName = strcat(inputFile,'_k_',num2str(K),'_n_',num2str(N),'_tau_',num2str(TAU));
%plotResults(fileName,data,bestC,bestMu,genome);

%Save Results to file
saveResults(fileName, bestC, bestMu, bestP, allLike);

resultsFile = strcat(fileName,'.results');

numResults = size(bestC,1);

plotResultsFromFile(resultsFile,inputFile,genome,11,numResults*4);

% Exit if flag is set - used for compiled matlab code
if (exitVal)
    exit;
end

end % end main function call


% This function parses the input parameters
function[N,K,TAU,exitVal] = parseInput(N,K,TAU,exitVal)

% Convert parameters to numbers (for running compiled code)
if (ischar(N))
    N=str2double(N);
end

if (ischar(K))
    K=str2double(K);
end

if (ischar(TAU))
    TAU=str2double(TAU);
end

if (ischar(exitVal))
    exitVal=str2double(exitVal);
end

end % end function call


% Get all C without regard to columns of all 0 or order of columns
% This function is now depreciated
function y = next(C,k)
    y = C;    
    for i=1:(size(C, 1)*size(C, 2))
        if(y(i) < k) 
           y(i) = y(i) + 1;
           for j=1:i-1
               y(j) = 0;                
           end
           return;
        end
    end
    y = zeros(size(C, 1), size(C, 2));
end % end function call


% Return the point on the simplex that corresponds to
% \widehat(WCmu)
function [P] = getSimplexPt(C,final_mu,normal)
    M = size(C,1);
    W1 = repmat(normal, 1, M);
    W2 = eye(M);
    W = W2 .* W1;

    WC = W*C;
    WCmu = WC*final_mu;

    P = WCmu./(sum(WCmu));
end

    
