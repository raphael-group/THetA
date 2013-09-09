function [hasNext, C] = generateNextC(aC, k, bounds, lowerBounds)

    hasNext = 1;
    C = aC;
    M = size(C, 1);
    N = size(C, 2);
    
    largest_col  = bounds;          

    for i=1:N
        if(lessThanCol(C(:, i), largest_col, M))
            %Now increase the i'th column
            for j=1:M
                if(C(j, i) < bounds(j))
                    C(j, i) = C(j, i) + 1;                
                    for h=1:(j-1)
                        C(h, i) = lowerBounds(h);
                    end  
                    break;
                end            
            end
            %making the previous column to smallest column
            for h=1:(i-1)
                C(:, h) = C(:, i);
            end
            return
        end        
    end
    %now here we figure out this was the last C. So put the first col 0        
    hasNext = 0;
end
