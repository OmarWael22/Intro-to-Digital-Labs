function sample_seq = GenerateSamples(bit_seq,fs)
%
% Inputs:
%   bit_seq:    Input bit sequence
%   fs:         Number of samples per bit
% Outputs:
%   sample_seq: The resultant sequence of samples
%
% This function takes a sequence of bits and generates a sequence of
% samples as per the input number of samples per bit

% repeat bit_seq in many rows
sample_seq = repmat(bit_seq',fs,1);
% combine rows in one row
sample_seq=sample_seq(:)';
% invert sample_seq to be one culomn
sample_seq=sample_seq';
