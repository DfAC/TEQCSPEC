function A = replace(A, S1, S2) ;
% REPLACE - Replace Elements
%   B = REPLACE(A,S1,S2) returns a matrix B in which the elements in A that 
%   are in S1 are replaced by those in S2. In general, S1 and S2 should have
%   an equal number of elements. If S2 has one element, it is expanded to
%   match the size of S1. Examples:
%      replace([1 1 2 3 4 4],[1 3],[0 99]) % ->  [ 0 0 2 99 4 4]
%      replace(1:10,[3 5 6 8],NaN) % ->  [ 1 2 NaN 4 NaN NaN 7 NaN 9 10]
%      replace([1 NaN Inf 8 99],[NaN Inf 99],[12 13 14]) % -> [1 12 13 8 14]
%
%   A and S1 can be cell arrays of strings. In that case S2 should be a
%   cell array as well but can contain mixed types. Example:
%      replace({'aa' 'b' 'c' 'a'},{'a' 'b'}, {'xx' 2}) %-> {'aa' [2] 'c' 'xx'}
%
%   If S2 is empty, the elements of A that are in S1
%   are removed. Examples:
%      replace(1:5,[2 4],[]) % -> [1 3 5]
%      replace({'aa' 'a' 'b' 'aa' 'c'},{'aa','c'},{}) % -> {'a', 'b'}
%
%   See also FIND, STRREP, REGEXPREP, ISMEMBER

% for Matlab R13
% version 1.2 (feb 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.0 (feb 2006) created
% 1.1 (feb 2006) fixed bug when NaNs were to be removed 
% 1.2 (feb 2006) fixed again bug with NaNs

error(nargchk(3,3,nargin)) ;

% all three inputs should be cell arrays or numerical arrays
if ~isequal(iscell(A), iscell(S1), iscell(S2)),
    error('The arguments should be all cell arrays or not.') ;
end

if iscell(A),
    % if they are cell, they should be character arrays
    if ~all(cellfun('isclass',A(:),'char')),
        error('A should be a cell array of strings.') ;
    end
    if ~all(cellfun('isclass',S1(:),'char')),
        error('S1 should be a cell array of strings.') ;
    end
end

if ~isempty(S2),
    if numel(S2)==1,
        % single element expansion
        S2 = repmat(S2,size(S1)) ;
    elseif numel(S1) ~= numel(S2),
        error('The number of elements in S1 and S2 do not match ') ;
    end
end

% the engine
[tf,loc] = ismember(A,S1) ;
if any(tf),
    if isempty(S2),
        A(tf) = [] ;
    else
        A(tf) = S2(loc(tf)) ;
    end
end

% special treatment for nans if necessary
if ~iscell(S1),
    % only for non-cell arrays
    qsn = isnan(S1(:)) ;
    if any(qsn),
        qa = isnan(A(:)) ;
        if any(qa),            
            if isempty(S2),
                A(qa) = [] ;
            else
                i = min(find(qsn)) ;            
                A(qa) = S2(i) ;
            end
        end
    end
end
