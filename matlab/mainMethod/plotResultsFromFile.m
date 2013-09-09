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
% This function plots the output data from the THetA program for mixtures
% of two or three genomes (n=2 or n=3).
%
% INPUT:
% resultsFile - the results file output by the THetA code.
% intervalFile - the input file require by THetA.
% genome - a name for the genome under consideration.
% width - the width of the pdf to create.
% height - the height of the pdf to create.
% 
% OPTIONAL PARAMETERS:
% rowNumber - plots only a particular entry (TO IMPLEMENT)
%
% OUTPUT:
% resultsFile.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = plotResultsFromFile(resultsFile, intervalFile, genome, width, height)

% Step 1: add access to main code
file_path = which('plotResultsFromFile.m');
f = filesep;
sep_locs = regexp(file_path,f);
idx = sep_locs(end - 1); %Go up one level
main = strcat(file_path(1:idx),'matlab/mainMethod');
addpath(genpath(main));

% WHERE TO PUT THIS CODE????

close all;
figure1 = figure;
hold on;

%Step 1: load interval file information

intervals = importdata(intervalFile);
%if (isstruct(intervals))
%    intervalMat = intervals.data;
%else
%    intervalMat = intervals;
%end

%If header line was included, we need to get the correct field
if (isfield(intervals,'data'))
    intervals = intervals.data;
end

%Check to see if ID column is still included

%if 6 columns then, remove the first one
if (size(intervals,2) == 6)
    intervals = intervals(:,2:end);
end

%if 8 columns then has bounds and ID, remove first
if (size(intervals,2) == 8)
    intervals = intervals(:,2:end);
end

intervalMat = intervals;

allData = importdata(resultsFile);
numResults = size(allData,1);

%Iterate over all results
for i=2:numResults  %start at 2 to ignore header file
    parts = regexp(allData{i}, '\t','split');
    likelihood = str2mat(parts(1));
    mu_str = regexp(parts{2},',','split');
    mu  = zeros(size(mu_str));
    mu = mu';
    for j=1:size(mu,1)
        mu(j) = str2num(mu_str{j});
    end
    
    
    c = [];
    c_str = regexp(parts{3},':','split'); 
    
    for j=1:size(c_str,2)
        c_row = regexp(c_str{j},',','split');
        row = zeros(1,size(c_row,2));
        for k=1:size(c_row,2)
            row(k) = str2num(c_row{k});
        end
        c = [c;row];
    end
        
        
    subplot(numResults - 1,1,i-1);
    hold on;
    plotTumor(intervalMat,mu,c, genome);
    
end

figuresize(width,height,'inches');
pdfName = strcat(resultsFile,'.pdf');
print(figure1,'-dpdf',pdfName);

end

function[] = plotTumor(intervalMat,mu,c,genome);
numIntervals = size(c,1);

legendPlot = zeros(size(c,2)+1,1);


vertBar = [-0.5,max(max(c))+ 0.5];
%offset = 0.075;
offset = (0.075 * max(max(c))/3);

tau= 2;

colors = {'b-','r-','k-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check for more components
if (size(c,2) > 2)
    fprintf('More than 3 components not currently supported');
    quit;
end

%Plot Tumor Bars
for i=1:size(c,2)
    pointer = 1;
    lastChrm = 1;
    lastEnd = 0;
    ticks = [1];
    chrms = [1];

    for k=1:numIntervals
        curChrm = intervalMat(k,1);
        curStart = intervalMat(k,2);
        curEnd = intervalMat(k,3);
        curLength = curEnd - curStart + 1;

        %Get chrm change and plot vert bar
        if (curChrm ~= lastChrm)
            marker = [pointer, pointer];
            plot(marker,vertBar,'k');
            ticks = [ticks,pointer];
            chrms = [chrms,curChrm];
        end

        %Check for skipped arm
        if (lastChrm ~= curChrm && curStart ~= 1)
            pointer = pointer + curStart;
        end
        
        %Check for skipped interval
        if (lastChrm == curChrm && lastEnd ~= curStart-1)
            pointer = pointer + (curStart - lastEnd);
        end

        xgrid = [pointer, pointer + curLength - 1];
        pointer = pointer + curLength;
        ygrid = [c(k,i)+ (i-1)*offset, c(k,i) + (i-1)*offset];
        p = plot(xgrid,ygrid,colors{i},'LineWidth',3);
        legendPlot(i) = p;

        lastChrm = curChrm;
        lastEnd = curEnd;

    end
    ticks = [ticks,pointer];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Normal Bars
pointer = 1;
lastChrm = 1;
lastEnd = 0;
ticks = [1];
chrms = [1];

for k=1:numIntervals
    curChrm = intervalMat(k,1);
    curStart = intervalMat(k,2);
    curEnd = intervalMat(k,3);
    curLength = curEnd - curStart + 1;

    %Get chrm change and plot vert bar
    if (curChrm ~= lastChrm)
        marker = [pointer, pointer];
        plot(marker,vertBar,'k');
        ticks = [ticks,pointer];
        chrms = [chrms,curChrm];
    end

    %Check for skipped arm
    if (lastChrm ~= curChrm && curStart ~= 1)
        pointer = pointer + curStart;
    end
    
    %Check for skipped interval
    if (lastChrm == curChrm && lastEnd ~= curStart-1)
        pointer = pointer + (curStart - lastEnd);
    end

    xgrid = [pointer, pointer + curLength - 1];
    pointer = pointer + curLength;
    %ygrid = [c(k,i)+ (i-1)*offset, c(k,i) + (i-1)*offset];
    ygrid = [tau - offset, tau - offset];
    p = plot(xgrid,ygrid,colors{3},'LineWidth',3);
    
    legendPlot(end) = p;

    lastChrm = curChrm;
    lastEnd = curEnd;

end
ticks = [ticks,pointer];
    


ylim(vertBar);
xlim([0,pointer]);
xlabel('Chromosome','FontSize',14);
ylabel('Copy Number','FontSize',14);

%Plot Chromosomes by name
revisedTicks = [];
for i=1:size(ticks,2)-1
    s1 = ticks(i);
    s2 = ticks(i+1);
    val = (s1 + s2)/2;
    revisedTicks = [revisedTicks,val];
end

set(gca,'XTick',revisedTicks);
set(gca,'XTickLabel',chrms);

if (size(c,2) == 2)
    h = legend(legendPlot,'Tumor1','Tumor2','Normal','Location','EastOutside');
    set(h,'FontSize',14);
else
    if (size(c,2) == 1)
        h=legend(legendPlot,'Tumor1','Normal','Location','EastOutside');
        set(h,'FontSize',14);
    else
        fprintf('More than 3 components not currently supported.');
        quit;
    end
end


figTitle = strcat(genome,' - Normal: ',num2str(100*mu(1),3),'%');

for j=2:size(mu,1)
    figTitle = strcat(figTitle, ', Tumor',num2str(j-1),': ',num2str(100*mu(j),3),'%'); 
end
title(figTitle, 'FontSize',14);

end

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
% Copyright and licence information appended.

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


