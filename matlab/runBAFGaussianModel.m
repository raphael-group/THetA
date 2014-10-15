%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will compute the NLL of the Gaussian model of BAF for the
% provided results file and will plot a figure of the best solution with
% the expected BAF based on the estimates of C and mu for the data on top
% of the BAF data.
%
% INPUT:
% 1. tumorSNP - tab delimited with the following columns for a set of known
%    germline SNP positions for the tumor BAM
%       1. Chrm
%       2. Pos
%       3. A - # reads with A at this position
%       4. C - # reads with C at this position
%       5. G - # reads with G at this position
%       6. T - # reads with T at this position
%       7. Total - # reads at this position
%       8. refCount - #reads with the reference allele
%       9. mutCount - #reads with the mutant allele
%
% 2. normalSNP - same format as tumorSNP but for the matched normal BAM
% 3. intervalFile -  The same interval file used when running THetA
% 4. resultsFile - a results file output by THetA
% 5. chrmsToUse (i.e. 1:22) - the list of chromosomes to consider (must be
%    numeric and match the format used in the other files
% 6. prefix - output prefix for created files
% 7. size - size of figure in inches (i.e. [11,8])
% 8. varargin - By default only a plot for the best C, mu pair is output.
%    Include an additional parameter 'all' to output plots of all solutions. 
%
% OUTPUT:
% prefix.BAF.BEST.pdf - pdf plot of best solution or prefix.BAF.ALL.pdf if
%                       the extra 'all' parameter is supplied
% prefix.BAF.NLL.results - a copy of the original results file with NLL
%                          appended for all solutions as an extra column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[NLL]=runBAFGaussianModel(tumorSNP, normalSNP, intervalFile, resultsFile, chrmsToUse, prefix, size, varargin)

%Step 0: Initialize any parameters or fixed variables in opts struct
opts.MIN_SIZE=2000000;
opts.MIN_SNP=10;
opts.MIN_HET=0.4;
opts.MAX_HET=0.6;
opts.TUM_READ=1;
opts.TYPE='best';

if nargin > 7
    opts.TYPE=lower(varargin{1});
end

%Step 1: Load in data
[tumor, normal] = loadSNPCountData(tumorSNP, normalSNP);
intervals = loadIntervalFile(intervalFile);
results = loadResultsFile(resultsFile);

%Step 2: calculate NLL for each
NLL = [];
for i=1:numel(results)
   C = results{i}.c;
   mu = results{i}.mu;
   
   curNLL = getNLL(tumor, normal, intervals, chrmsToUse,C,mu,opts);
   NLL = [NLL;curNLL];
   results{i}.BAF_NLL=curNLL;
end

%Step 3: Save NLL to file
saveResultsWithNLL(resultsFile,prefix,NLL);

%Step 4: Create plot
plotBAFGaussianResults(tumor, normal, intervals, chrmsToUse, results, prefix,size,opts);


end

% Plot all the BAFs points for heterozygous SNPs
function[h1]=plotBAF(tumor, normal, chrmsToUse, opts)

MIN_HET=opts.MIN_HET;
MAX_HET=opts.MAX_HET;
TUM_READ=opts.TUM_READ;

tumorBAF = tumor(:,9)./(tumor(:,8) + tumor(:,9));
normalBAF = normal(:,9)./(normal(:,8) + normal(:,9));


minR = 0;
maxR = 1;

vertBar = [minR, maxR];

offset = 0;
ticks = [0];
barX = [];
barY = [];

chrms = tumor(:,1);
uniqueVal = unique(chrms);
numChrms = numel(uniqueVal);
list = hsv(numChrms);

legendVal = cell(size(chrmsToUse,2),1);

hetSNP = find(normalBAF > MIN_HET & normalBAF < MAX_HET);
posReads = find(tumor(:,8) >= TUM_READ | tumor(:,9) >= TUM_READ);

for i=1:numel(uniqueVal)
   chrm = uniqueVal(i);
   
   if (~isempty(find(chrmsToUse == chrm)))
       rows_temp = find(chrms == chrm);
       temp1 = intersect(rows_temp, hetSNP);
       rows = intersect(temp1, posReads);

       xCoords = tumor(rows,2) + offset;
       
       %Update offset by chrm
       all_xCoords = tumor(rows_temp,2);
       offset = offset + max(all_xCoords);
       
       %offset = max(xCoords);
       ticks = [ticks,offset];
       yCoords = tumorBAF(rows);

       hold on;
       color = list(i,:);
       h1 = plot(xCoords,yCoords,'.','Color',color);

       idx = find(chrmsToUse == chrm);

       legendVal{idx} = strcat('Chr',num2str(chrm));

       barX = [barX;offset,offset];
       barY = [barY;vertBar];
   end
end

plot(barX',barY','k','LineWidth',1);
ylim(vertBar);
xlim([0,offset]);


xlabel('Chromosome','FontSize',14);
ylabel('BAF','FontSize',14);

%lko - 3/16/2014 - comment out legend
%h = legend(legendVal);
%set(h,'FontSize',14);
%a=get(h,'children');

%a(1) corresponds to the marker object
%set(a([1:3:end]),'markersize',25);

%Plot Chromosome by name
revisedTicks = [];
for i=1:size(ticks,2)-1
    s1=ticks(i);
    s2=ticks(i+1);
    val = (s1 + s2)/2;
    revisedTicks = [revisedTicks,val];
end

set(gca,'XTick',revisedTicks);
set(gca,'XTickLabel',chrmsToUse);
    
end


function[h1] = plotOurEstimate(tumor, normal, intervals, chrmsToUse,C,mu,opts)

%Fixed parameters
MIN_SIZE=opts.MIN_SIZE;
MIN_SNP=opts.MIN_SNP;
MIN_HET=opts.MIN_HET;
MAX_HET=opts.MAX_HET;
TUM_READ=opts.TUM_READ;

tumorBAF = tumor(:,9)./(tumor(:,8) + tumor(:,9));
normalBAF = normal(:,9)./(normal(:,8) + normal(:,9));

hetSNP = find(normalBAF > MIN_HET & normalBAF < MAX_HET);
posReads = find(tumor(:,8) >= TUM_READ | tumor(:,9) >= TUM_READ);

numChrms = numel(chrmsToUse);
offset = 0;

for i=1:numChrms
    chrm = chrmsToUse(i);
  
    %Go through each interval
    int_rows = find(intervals(:,1) == chrm);
       
   for j=1:numel(int_rows)
       start_pos = intervals(int_rows(j),2);
       end_pos = intervals(int_rows(j),3);
       length = end_pos - start_pos +1;

       %Get all BAF in range
       chrm_rows = find(tumor(:,1) == chrm);
       start_rows = find(tumor(:,2) >= start_pos);
       end_rows = find(tumor(:,2) <= end_pos);
       temp1 = intersect(chrm_rows, start_rows);
       temp2 = intersect(hetSNP,temp1);
       temp3 = intersect(temp2, posReads);
       rows = intersect(temp3,end_rows);

       if (numel(rows > MIN_SNP) && (length > MIN_SIZE))
            xCoords = tumor(rows,2);
            tumor_yCoords = tumorBAF(rows);
            %Get rid of any NaN
            nan_idx = find(~isnan(tumor_yCoords));
            tumor_yCoords = tumor_yCoords(nan_idx);
            
            normal_yCoords = normalBAF(rows);
            
            tumor_std = std(tumor_yCoords);
            normal_std = std(normal_yCoords);
            
            start_plot = start_pos + offset;
            end_plot = offset + start_pos + length;

            xVals = [start_plot, end_plot];
            
            [above, below] = calcEstBAF(C,mu,int_rows(j));
            
            if (above == below)
                h1 = plot(xVals, [above, above],'-k','LineWidth',2);
            else
                h1 = plot(xVals, [above, above],'-k','LineWidth',2);
                h1 = plot(xVals, [below, below],'-k','LineWidth',2);
            end
            
       end
   end
   
   %Update offset by chrm
   BAF_chrm = find(tumor(:,1) == chrm);
   xCoords = tumor(BAF_chrm,2);
   offset = offset + max(xCoords);
end
   


end


function [above, below] = calcEstBAF(C,mu,idx)

n = size(C,2);

%invalid = find(C(idx,:) < 1 | C(idx,:) > 3); %Allow 0 copies
invalid = find(C(idx,:) > 3);

%If copy 0 or above 3 we can't estimate parental alleles.
if (numel(invalid) > 0)
    above = NaN;
    below = NaN;
    return;
end

denom = sum(C(idx,:)*mu);
above_num = 0;
below_num = 0;

for i=1:n
   
    if (C(idx,i) == 2)
        above_num = above_num + mu(i);
        below_num = below_num + mu(i);
    end
    
    if (C(idx,i) == 1)
        above_num = above_num + mu(i);
    end
    
    if (C(idx,i) == 3)
        above_num = above_num + 2*mu(i);
        below_num = below_num + mu(i);
    end
    
end

above = above_num/denom;
below = below_num/denom;

end

function [NLL] = getNLL(tumor, normal, intervals, chrmsToUse,C,mu,opts)

%FIXED PARAMETERS
MIN_SIZE=opts.MIN_SIZE;
MIN_SNP=opts.MIN_SNP;
MIN_HET=opts.MIN_HET;
MAX_HET=opts.MAX_HET;
TUM_READ=opts.TUM_READ;

NLL=0;

tumorBAF = tumor(:,9)./(tumor(:,8) + tumor(:,9));
normalBAF = normal(:,9)./(normal(:,8) + normal(:,9));

hetSNP=find(normalBAF >= MIN_HET & normalBAF <= MAX_HET);

numChrms = numel(chrmsToUse);

for i=1:numChrms
    chrm = chrmsToUse(i);
  
    %Go through each interval
    int_rows = find(intervals(:,1) == chrm);
       
   for j=1:numel(int_rows)
       start_pos = intervals(int_rows(j),2);
       end_pos = intervals(int_rows(j),3);
       length = end_pos - start_pos +1;

       %Get all BAF in range
       chrm_rows = find(tumor(:,1) == chrm);
       start_rows = find(tumor(:,2) >= start_pos);
       end_rows = find(tumor(:,2) <= end_pos);
       pos_reads = find(tumor(:,8) >= TUM_READ | tumor(:,9) >= TUM_READ);
       temp1 = intersect(chrm_rows, start_rows);
       temp2 = intersect(hetSNP,temp1);
       temp3 = intersect(temp2, pos_reads);
       rows = intersect(temp3,end_rows);

       if (numel(rows > MIN_SNP) && (length > MIN_SIZE))
            xCoords = tumor(rows,2);
            tumor_yCoords = tumorBAF(rows);
            
            %Get rid of any NaN for exome data (already handled now by at
            %least one read
            nan_idx = find(~isnan(tumor_yCoords));
            tumor_yCoords = tumor_yCoords(nan_idx);
            
            normal_yCoords = normalBAF(rows);
            
            %tumor_std = std(tumor_yCoords); %not used - comment out
            normal_std = std(normal_yCoords);
            
            above_y_idx = find(tumor_yCoords >= 0.5);
            above_y = tumor_yCoords(above_y_idx);
            below_y_idx = find(tumor_yCoords < 0.5);
            below_y = tumor_yCoords(below_y_idx);
            
            [above, below] = calcEstBAF(C,mu,j);
            above_mu_vector = repmat(above,numel(above_y),1);
            below_mu_vector = repmat(below,numel(below_y),1);
            std_above_vector = repmat(normal_std,numel(above_y),1);
            std_below_vector = repmat(normal_std, numel(below_y),1);
            
            
            if (numel(above_y) > MIN_SNP)
                NLL = NLL + sum(-log(normpdf(above_y,above_mu_vector, std_above_vector)));
            end
            
            
            if (numel(below_y) > MIN_SNP)
                NLL = NLL + sum(-log(normpdf(below_y,below_mu_vector, std_below_vector)));
            end
            
       end%end check for long enough and enough SNPs
   end%end iterate on intervals
end %end check all chrms

end

%Loads in SNP file data from tumor and normal data.  Assumes that all data
%in the files is numeric.  Returns the data as matrices
function [tumor, normal] = loadSNPCountData(tumorSNP, normalSNP)

tumor = importdata(tumorSNP);

if (isstruct(tumor))
    tumor = tumor.data;
end

normal = importdata(normalSNP);

if (isstruct(normal))
    normal = normal.data;
end

end %end loadSNPCountData function


%Loads in data from an interval file.  Handles the fact that the first
%column (ID) may be numeric or not.  Does so by only looking at the number
%of columns loaded in.  Beware if the interval file also has lower and
%upper bound information this may be problematic!
function[intervals] = loadIntervalFile(intervalFile)

intervals = importdata(intervalFile);

if (isstruct(intervals))
    intervals = intervals.data;
    
    if (size(intervals,2) == 6) %ID column still there
        intervals = intervals(:,2:end);
    end
end

end %end loadIntervalFile function


%Loads in the data from a THetA results file into a list of structures with
%fields C and mu.
function[results] = loadResultsFile(resultsFile)

    %fixed params
    tau = 2;

    %new way to import data
    fid = fopen(resultsFile);
    header = fgetl(fid);
    allData = {header};
    line = fgetl(fid);

    while (line ~= -1)
        allData{numel(allData)+1} = line;
        line = fgetl(fid);
    end
    fclose(fid);

    numResultsRows = size(allData,2);
    results = cell(numResultsRows - 1,1);

    %Iterate over all results
    for i=2:numResultsRows  %start at 2 to ignore header file

        parts = regexp(allData{i}, '\t','split');

        if (numel(parts) == 1)
            continue;
        end

        likelihood = str2mat(parts(1));
        mu_str = regexp(parts{2},',','split');
        mu  = zeros(size(mu_str));
        mu = mu';
        for j=1:size(mu,1)
            mu(j) = str2num(mu_str{j});
        end


        c = [];
        c_str = regexp(parts{3},':','split'); 
        neg_idx = [];

        for j=1:size(c_str,2)
            c_row = regexp(c_str{j},',','split');
            row = zeros(1,size(c_row,2));
            for k=1:size(c_row,2)

                if (strcmp(c_row{k},'X'))
                    row(k) = -1;
                    neg_idx = [neg_idx;j];
                else
                    row(k) = str2num(c_row{k});
                end
            end
            c = [c;row];
        end

        %Check if we should flip columns to make components in reverse sorted order
        if (numel(mu) == 3 && mu(2) < mu(3))
            temp = mu(3);
            mu(3) = mu(2);
            mu(2) = temp;
            temp = c(:,1);
            c(:,1) = c(:,2);
            c(:,2) = temp;
        end

        %Add normal to c
        num = size(c,1);
        normal = tau * ones(num,1);

        %Change any normal in row with -1 (not estimated to -1)
        normal(neg_idx) = -1;

        c = [normal, c];
        
        obj.c = c;
        obj.mu = mu;
        results{i-1} = obj;
    end
end %end loadResultsFile function


%This function will save a new results file with the computed NLL for the
%BAF gaussian generative model as a new column.  The file is save with the
%specified prefix.
function [] = saveResultsWithNLL(resultsFile,prefix,NLL)

    %new way to import data
    fid = fopen(resultsFile);
    header = fgetl(fid);
    allData = {header};
    line = fgetl(fid);

    while (line ~= -1)
        allData{numel(allData)+1} = line;
        line = fgetl(fid);
    end
    fclose(fid);
    
    %Append new header
    allData{1} = sprintf('%s\t%s',allData{1},'BAF_NLL');
    
    %Append NLL stuff
    for i=2:numel(allData)
       allData{i} = sprintf('%s\t%.2f',allData{i},NLL(i-1)); 
    end
    
    %Save to file
    file=strcat(prefix,'.BAF.NLL.results');
    fid=fopen(file,'w');
    
    for i=1:numel(allData)
       fprintf(fid,'%s\n',allData{i}); 
    end
    fclose(fid);
end

%This function will call the plotting code for the program.
function[] = plotBAFGaussianResults(tumor, normal, intervals, chrmsToUse, results, prefix,size,opts)

%Step 0:setup plotting and get type
close all;
figure1=figure;
hold on;
fileName = strcat(prefix,'.BAF.BEST.pdf');

%Step 1: go through and do plotting
if strcmp(opts.TYPE,'all')
    numPlot = numel(results);
    for i=1:numPlot
        subplot(numPlot,1,i);
        C=results{i}.c;
        mu=results{i}.mu;
        NLL=results{i}.BAF_NLL;
        plotWholeResult(tumor, normal, chrmsToUse, intervals, C, mu, NLL,opts)
    end
    fileName = strcat(prefix,'.BAF.ALL.pdf');
else
    %Find best solution
    numPlot = numel(results);
    bestNLL = inf;
    bestIdx = 1;
    for i=1:numPlot
        NLL=results{i}.BAF_NLL;
        if NLL < bestNLL
            bestIdx = i;
            bestNLL = NLL;
        end
    end
    
    C=results{bestIdx}.c;
    mu=results{bestIdx}.mu;
    NLL=results{bestIdx}.BAF_NLL;
    plotWholeResult(tumor, normal, chrmsToUse, intervals, C, mu, NLL,opts)
    
end %end if type


%Step 2: save the plot
figuresize(size(1),size(2),'inches');
print(figure1,'-dpdf','-r150',fileName);



end


%Pull all the plotting stuff into one method
function [] = plotWholeResult(tumor, normal, chrmsToUse, intervals, C, mu, NLL,opts)
    h1 = plotBAF(tumor, normal, chrmsToUse,opts);
    h2 = plotOurEstimate(tumor, normal, intervals, chrmsToUse,C,mu,opts);
    
    titleStr=sprintf('BAF Model NLL: %.2f',NLL);
    title(titleStr);
    
    legend_handle = legend([h1 h2],'Observed BAF','Expected BAF');
    
    %Force legend colors to be black (since legend is just for shape)
    leg_line=findobj(legend_handle,'type','Line');
    for i = 1:length(leg_line)
        set(leg_line(i), 'Color', 'k');
    end
    
end

function figuresize( w , h , u )
%FIGURESIZE Set a figure to a specific size
%
% When saving a figure as a PDF, it is necessary to set the
% figure size appropriately. This function sets the "paper size"
% sufficient that the figure is saved with a tight bounding box.
% It will also set the figure size on screen correspondingly.
%
% figuresize(width,height) - sets the figure size in centimeters
% figuresize(width,height,units) - sets the figure size in <units>
%
% <units> can be any of the standard Matlab lengths.
%
% Will Robertson
% 28 April 2010
%
% Copyright and licence information appended below.
%
% Copyright (c) 2010 Will Robertson, wspr 81 at gmail dot com
% All rights reserved.
%
% Distributed under the BSD licence in accordance with the wishes of the
% Matlab File Exchange. (Usually I'd pick the Apache License.)
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ''AS IS'' AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

p = inputParser;
p.addRequired('width', @(x) isnumeric(x) && all(size(x)==1) );
p.addRequired('height',@(x) isnumeric(x) && all(size(x)==1) );
p.addOptional('units','centimeters',...
  @(x) any(strcmp(x,{'normalized','centimeters','inches','points'})) );

p.parse( w, h, u );
w = p.Results.width;
h = p.Results.height;
u = p.Results.units;

p = 0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
  'Position',[screenpos(1:2) w h],...
  'PaperUnits',u,...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);

end


