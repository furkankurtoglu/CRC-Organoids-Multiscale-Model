
function [h, adj_p]=fdr_bonferroni(pvals,q,n_tests)
    if nargin<1,
        error('You need to provide a vector or matrix of p-values.');
    else
        if ~isempty(find(pvals<0,1)),
            error('Some p-values are less than 0.');
        elseif ~isempty(find(pvals>1,1)),
            error('Some p-values are greater than 1.');
        end
    end
    if nargin<2,
        q=.05;
    end
    if nargin<3,
        n_tests = numel(pvals);
    end

    adj_p = pvals * n_tests;
    adj_p = min(adj_p, 1);
    h = (adj_p <= q);
end
