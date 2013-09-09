function [fbC, fbMu, fbP, fbLike] = filterBests(bestC, bestMu, bestP, allLike, bestLikelihood, epsilon)
    
    indicesToRemove = [];
    
    fbC    = bestC;
    fbMu   = bestMu;
    fbP    = bestP;
    fbLike = allLike;
    
    for i=1:length(allLike)
        if(allLike(i) > bestLikelihood + epsilon)
            indicesToRemove = [indicesToRemove, i - length(indicesToRemove)];
        end
    end
    
    for i = indicesToRemove
        fbC(i) = [];
        fbMu(i) = [];
        fbP(i) = [];
        fbLike(i) = [];
    end
    
end
