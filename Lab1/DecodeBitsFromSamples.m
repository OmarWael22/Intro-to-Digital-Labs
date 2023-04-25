function rec_bit_seq = DecodeBitsFromSamples(rec_sample_seq,case_type,p,fs)
%
% Inputs:
%   rec_sample_seq: The input sample sequence to the channel
%   case_type:      The sampling frequency used to generate the sample sequence
%   fs:             The bit flipping probability
% Outputs:
%   rec_sample_seq: The sequence of sample sequence after passing through the channel
%
% This function takes the sample sequence after passing through the
% channel, and decodes from it the sequence of bits based on the considered
% case and the sampling frequence

if (nargin <= 2)
    fs = 1;
end

switch case_type
    
    case 'part_1'
        % creating decide vector
        Decide = rand(size(rec_sample_seq))<=p;
        % calculate sequence after receiver
        rec_bit_seq = xor(rec_sample_seq, Decide);
    case 'part_2'
        % creating decide vector
        Decide = rand(size(rec_sample_seq))<=p;
        % calculate sequence after receiver before merging fs values
        rec_bit = xor(rec_sample_seq, Decide);
        rec_bit_seq = zeros(length(rec_sample_seq)/fs,1);
        for i=0:(length(rec_bit_seq)-1)
            % sum of each 'fs' bits
            sum=0;
            for j=1:fs
                sum = sum + rec_bit(i*fs+j);
            end
            % check if 1's is repeated more than the half of the 'fs' bits
            if(sum>=(fs/2))
                rec_bit_seq(i+1) = 1;
            end
        end
            
    case 'part_3'
        % creating decide vector
        Decide = rand(size(rec_sample_seq))<=p;
        % calculate sequence after receiver before merging fs values
        rec_bit = xor(rec_sample_seq, Decide);
        rec_bit_seq = zeros(length(rec_sample_seq)/fs,1);
        for i=0:(length(rec_bit_seq)-1)
            % sum of each 10 bits
            sum=0;
            for j=1:10
                sum = sum + rec_bit(i*fs+j);
            end
            % check if 1's is repeated more than the half of the 10 bits
            if(sum>=(fs/2))
                rec_bit_seq(i+1) = 1;
            end
        end
end

        
        
        
        