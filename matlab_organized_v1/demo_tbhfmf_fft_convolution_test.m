% Compare FFT convolution with correct direct cyclic sum
% for CENTERED-order storage (q=0 at center).

N = 17;
rng(1);

ix = (0:N-1) - floor(N/2);
[iy,ix2] = ndgrid(ix,ix);
qabs = sqrt(ix2.^2 + iy.^2) + 1e-17;

Vq = 1 ./ qabs;
c  = floor(N/2)+1;
% Vq(c,c) = 10;

rho = randn(N,N) + 1i*randn(N,N);

% ----- FFT convolution (centered -> FFT order -> fft2) -----
to_fft   = @(A) ifftshift(ifftshift(A,1),2);
from_fft = @(A) fftshift(fftshift(A,1),2);

Vr   = fft2(to_fft(Vq));
rhoR = fft2(to_fft(rho));
Sigma_fft = from_fft(ifft2(Vr .* rhoR));  % centered order output

% ----- Direct cyclic sum (CENTERED kernel indexing) -----
Sigma_dir = zeros(N,N);
for iky = 1:N
  for ikx = 1:N
    s = 0;
    for iqy = 1:N
      for iqx = 1:N
        % IMPORTANT: iq corresponds to shift (iq-c), not (iq-1)
        jky = mod( (iky - iqy + c - 1), N ) + 1;
        jkx = mod( (ikx - iqx + c - 1), N ) + 1;
        s = s + Vq(iqy,iqx) * rho(jky,jkx);
      end
    end
    Sigma_dir(iky,ikx) = s;
  end
end

err = max(abs(Sigma_fft(:) - Sigma_dir(:)));
fprintf('max|Sigma_fft - Sigma_dir| = %.3e\n', err);

figure; imagesc(real(Sigma_fft)); axis image; title('Re Sigma (FFT)'); colorbar;
figure; imagesc(real(Sigma_dir)); axis image; title('Re Sigma (direct, fixed)'); colorbar;

