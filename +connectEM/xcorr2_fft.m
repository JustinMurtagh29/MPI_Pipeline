function [x_shift,y_shift] = xcorr2_fft(A, B)
    % Normalized FFT-based cross-correlation for N dimensional arrays

    % NOTE: This is not normalized cross correlation as defined on wikipedia:
    % https://en.wikipedia.org/wiki/Cross-correlation#Normalized_cross-correlation
    % Rather each image is normalized to zero mean and unit STD, maybe that is what is meant in the paper
    A = single(A);
    B = single(B);
    A = (A - mean(A(:))) ./ std(A(:));
    A = (B - mean(B(:))) ./ std(B(:));

    nd = max(ndims(A),ndims(B));
    dims = 1:nd;
    dims = reshape(dims, 1, []);

    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        l = m+n-1;
        A = fft(A,l,dim);
        B = fft(B,l,dim);
        subs{dim} = 1:m+n-1;
    end
    A = A.*conj(B);
    for dim=dims
        A = ifft(A,[],dim);
    end
    A = real(A(subs{:}));
    [x_max, y_max] = find(A == max(A(:)));
    x_shift = x_max - size(A,1) -1;
    y_shift = y_max - size(A,2) -1;

end

