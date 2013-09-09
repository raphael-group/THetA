%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a projection of the supplied mu onto the 
% probability simplex.
% 
% INPUT: 
% mu - found vector of mu.
%
% OUTPUT:
% proj_mu -  The mu values are mapped to the closest point on the simplex 
% using the L2 distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[proj_mu] = projectSimplex(mu)

N = size(mu,1);

% Project mu onto the simplex
temp = sort(mu,'descend');

sm = 0;
for j = 1:N
    sm = sm + temp(j);
    if temp(j) - (1/j)*(sm-1) > 0
        row = j;
        sm_row = sm;
    end
end
theta = (1/row)*(sm_row-1);
proj_mu = max(mu-theta,0);