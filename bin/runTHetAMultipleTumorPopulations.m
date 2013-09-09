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
function[bestLike, bestMu, bestC, bestP, allLike] = runTHetAMultipleTumorPopulations(inputFile, N, K, TAU, genome)

% Step 1: add access to main code
file_path = which('runTHetAMultipleTumorPopulations.m');
f = filesep;
sep_locs = regexp(file_path,f);
idx = sep_locs(end - 1); %Go up one level
main = strcat(file_path(1:idx),'matlab/mainMethod');
addpath(genpath(main));

runTHetA(inputFile, N, K, TAU, genome, 0);

end