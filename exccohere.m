function EC = exccohere(CptC,varargin)
%EXCCOHERE  Calculates excluded coherence (EC) for inputs or subsets of 
% inputs in a given multi-input, single output system
%   EC = EXCCOHERE(CptC), returns excluded coherence for each input,
%   calculated as a sum of component coherence terms. Component coherence
%   is calculated using cptcohere, see "help cptcohere" for complete
%   details. Excluded coherence describes the portion of the output power
%   that would be removed if the given input (or input subset) were
%   excluded. EC is a q x n array, where n is the length of the first
%   dimension of CptC, and q is the number of inputs to the system. EC(f,i)
%   is the excluded coherence between input i and the output at frequency 
%   f. 
%
%   EC = EXCCOHERE(CptC,S), returns excluded coherence for the input 
%   subsets defined in S. Each column of S describes an input subset for 
%   which excluded coherence is to be calculated. EC is returned as an 
%   array with a separate column for excluded coherence between each input 
%   subset and the output, in the same order specified by S. If subsets are
%   of different sizes, pad columns in S with NaN for smaller subsets. 
%   Example: For a given 4-input system, it is of interest to calculate
%   excluded coherence for 5 input subsets: inputs 1,2,3; inputs 1,3,4; 
%   inputs 1,3; inputs 2,4; input 4. S for this case should be
%       S = [1    1    1    2    4;
%            2    3    3    4    NaN;
%            3    4    NaN  NaN  NaN];
%
% Written by Dr. Steven Charles and Nolan Howes, 2024.


    % Manage additional inputs
    narginchk(1,2)
    q = length(CptC(1,1,:));
    S = linspace(1,q,q);
    if nargin > 1; S = varargin{1}; end

    % Calculate excluded coherence by summing correct component coherence
    % terms
    EC = zeros(length(CptC(:,1,1)),length(S(1,:)));
    for i = 1:length(S(1,:)) %for each column in S
        for j = 1:q
            for k = 1:q
                if ismember(j,S(:,i)) || ismember(k,S(:,i)) %if j OR k is a member of the subset
                    EC(:,i) = EC(:,i) + squeeze(CptC(:,j,k)); %Add the given component coherence term
                end
            end
        end
    end
    EC = real(EC); %Removes zero-valued imaginary part
end
    




