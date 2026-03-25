function [MH, invMH, M, invM, low_idx_0] = build_low_freq_projection(N, H)

low_indx  = 0:(H-1)/2;

high_indx = N-(H-1)/2:N-1;

low_idx_0 = [low_indx, high_indx];
low_idx   = low_idx_0 + 1;

row_indices = reshape([2*low_idx-1; 2*low_idx],1,[]);

M = calc_M_real(N);
MH = M(row_indices,:);

invM = (1/N)*M';
invMH = (1/N)*MH';

end
