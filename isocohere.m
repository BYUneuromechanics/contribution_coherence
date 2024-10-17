function IC = isocohere(CptC,varargin)
%ISOCOHERE  Calculates isolated coherence (IC) for inputs or subsets of 
% inputs in a given multi-input, single output system
%   IC = ISOCOHERE(CptC), returns isolated coherence for each input,
%   calculated as a sum of component coherence terms. Component coherence
%   is calculated using cptcohere, see "help cptcohere" for complete
%   details. Isolated coherence describes the portion of the output power
%   that would remain if the given input (or input subset) were isolated.
%   IC is a q x n array, where n is the length of the first dimension of
%   CptC, and q is the number of inputs to the system. IC(f,i) is the
%   isolated coherence between input i and the output at frequency f.
%
%   IC = ISOCOHERE(CptC,S), returns isolated coherence for the input 
%   subsets defined in S. Each column of S describes an input subset for 
%   which isolated coherence is to be calculated. IC is returned as an 
%   array with a separate column for isolated coherence between each input 
%   subset and the output, in the same order specified by S. If subsets are 
%   of different sizes, pad columns in S with NaN for the smaller subsets. 
%   Example: For a given 4-input system, it is of interest to calculate
%   isolated coherence for 5 input subsets: inputs 1,2,3; inputs 1,2,4; 
%   inputs 1,2; inputs 1,3; input 1. S for this case should be
%       S = [1    1    1    1    1;
%            2    2    2    3    NaN;
%            3    4    NaN  NaN  NaN];
%
% Written by Dr. Steven Charles and Nolan Howes, 2024.

    % Manage additional inputs
    narginchk(1,2)
    q = length(CptC(1,1,:));
    S = linspace(1,q,q);
    if nargin > 1; S = varargin{1}; end

    % Calculate isolated coherence by summing correct component coherence
    % terms
    IC = zeros(length(CptC(:,1,1)),length(S(1,:)));
    for i = 1:length(S(1,:)) %for each column in S
        for j = 1:q
            for k = 1:q
                if ismember(j,S(:,i)) && ismember(k,S(:,i)) %if j AND k are members of the subset
                    IC(:,i) = IC(:,i) + squeeze(CptC(:,j,k)); %Add the given component coherence term
                end
            end
        end
    end
    IC = real(IC); %Removes zero-valued imaginary part
end
    




