%Checks to see if all pairs of rows have some elements in the correct
%sorted order as determined by the sorted_idx vector (increasing)
function[isValid] = checkValidC(C,sorted_idx)

M = size(C,1);
isValid=1;

for i=1:M
   idx1=sorted_idx(i);
   curRow = C(idx1,:);
   
   for j=1:i-1
       idx2=sorted_idx(j);
       otherRow = C(idx2,:);
       
       %Make sure some column has idx2 <= idx1
       comp = (otherRow <= curRow);
    
        idx = find(comp == 1);
    
        if (isempty(idx))
            isValid = 0;
            return;
        end
   end
end
