function BER = ComputeBER(bit_seq,rec_bit_seq)
%
% Inputs:
%   bit_seq:     The input bit sequence
%   rec_bit_seq: The output bit sequence
% Outputs:
%   BER:         Computed BER
%
% This function takes the input and output bit sequences and computes the
% BER
% number of error bits
sum = 0;
for i = 1:size(rec_bit_seq)
    if(bit_seq(i)~=rec_bit_seq(i))
        sum= sum + 1;
    end
end
%calculate bit error rate
BER = sum / length(bit_seq);
