%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the multinomial likelihood, gradient, and hessian
% for the supplied values.
%
% INPUT:
% mu
% C (full matrix, with normal column)
% tumor
% normal
% N
%
% OUTPUT:
% likelihood
% gradient
% hessian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[likelihood, gradient] = multinomialLikelihood(mu, C, tumor, normal, N)

M = size(tumor,1);

% Compute useful matrices
W1 = repmat(normal, 1, M);
W2 = eye(M);
W = W2 .* W1;

WC = W*C;
WCmu = WC*mu;

%likelihood = 0;
%for i=1:M
%   likelihood = likelihood - tumor(i)*log(WCmu(i)/sum(WCmu)); 
%end

%Vectorize likelihood
hatWCmu = WCmu./sum(WCmu);
loghatWCmu = log(hatWCmu);
likelihood = -tumor'*loghatWCmu;


if (nargout == 2)
% Do gradient
    gradient = zeros(N,1);
    for h=1:N

        %for i=1:M
        %   gradient(h) = gradient(h) + -tumor(i)*(WC(i,h)/WCmu(i)); 
        %   gradient(h) = gradient(h) + tumor(i)*(sum(WC(:,h)/sum(WCmu)));
        %end

        col1 = WC(:,h)./WCmu;
        p1 = -tumor' * col1;

        num = sum(WC(:,h));
        val = num/sum(WCmu);
        col2 = repmat(val,M,1);
        p2 = tumor'*col2;
        gradient(h) = p1 + p2;

    end
end

% Do Hessian
%hessian = zeros(N,N);
%for h=1:N
%    for g=1:N
%        
%        for i=1:M
%           hessian(h,g) = hessian(h,g) + tumor(i)* (WC(i,h)*WC(i,g)/(WCmu(i)^2));
%           hessian(h,g) = hessian(h,g) + -tumor(i) *(sum(WC(:,h))*sum(WC(:,g))/(sum(WCmu))^2);
%        end
%        
%    end
%end

