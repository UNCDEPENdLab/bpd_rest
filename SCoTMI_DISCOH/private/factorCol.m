% [selectedFactor, selectedCols] = 
%               factorCol(matrix, factorNumber, factorSize)
%
%   Chooses groups of columns from a matrix. Column groups are of uniform
%   size. Can be used for row vectors (or transposed column vectors) to
%   extract groups of elements.
% 
% Inputs:
%   factorSize: vector of sizes for each sub-group, in order
%   factorNumber: which sub-group we want
%   matrix: to choose from, size(matrix)==[whatever, sum(factorSize)]
%
% Outputs:
%   selectedFactor: 'factorNumber' column-group of matrix
%   selectedCols: selected columns
%
%   Ben Cassidy 22/2/2011
%

function [selectedFactor, selectedCols] = factorCol(matrix, factorNumber,factorSize)
    j=factorNumber;
    p=factorSize;
    selectedCols = (sum(p(1:(j-1)))+1 : sum(p(1:(j))));
    selectedFactor = matrix(:,selectedCols);
    
end
