function M = calc_M_real(N)
% Build the 2N x 2N real block DFT matrix.
%
% The complex DFT matrix entry W(k,n) = exp(-i*2*pi*k*n/N)
% is replaced by the real 2x2 block:
%
%   [ cos(-theta)   -sin(-theta)
%    sin(-theta)   cos(-theta) ]
%
% where theta = -2*pi*k*n / N.
%
% Output M is a real matrix of size (2N x 2N).

    M = zeros(2*N, 2*N);  % pre-allocate space

    for k = 0:N-1
        for n = 0:N-1
            
            theta = -2*pi*k*n / N;

            % 2x2 rotation block
            B = [ cos(theta),  -sin(theta);
                 sin(theta),  cos(theta) ];

            % Insert block into the right place
            rowIdx = 2*k + (1:2);
            colIdx = 2*n + (1:2);
            M(rowIdx, colIdx) = B;
        end
    end
end