function A = zscore(X)
    A = X - mean(X, 2);
    A = bsxfun(@rdivide, A, std(A, [], 2));    
end