function y = lessThanCol(a, b, M)
    y = 1;
    i = M;
    while(i >= 1) 
        if(a(i) == b(i))
            i = i-1;
        elseif (a(i) < b(i))
            y = 1;
            return;
        else
            y = 0;
            return;
        end
    end
    y=0;
end