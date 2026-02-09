function [Dk_core, gk] = get_Dk_qmd_plainD_core( ...
        Hk, Unk, Enk, band_list, dk_vecs, deltaE_reg)
% GET_DK_QMD_PLAIND_CORE_WITH_G
%   在 2D k-mesh 上（周期边界中心差分）计算并保存：
%
%   1) plain-QMD 核 (k-resolved):
%       D^{(n)}_{abc}(k) = v_a^n(k) g^{(n)}_{bc}(k) - v_b^n(k) g^{(n)}_{ac}(k)
%
%   2) quantum metric (k-resolved):
%       g^{(n)}_{ab}(k) = Re Σ_{m≠n} M^a_{nm}(k) M^b_{mn}(k) / (ε_n-ε_m)^2
%
%   M^a_{nm}(k) = <u_n(k)| ∂_{k_a} H(k) |u_m(k)>,  a=x,y
%   v_a^n(k)    = (1/ħ) Re M^a_{nn}(k)
%
% Inputs:
%   Hk        : nb x nb x Nkx x Nky, H(k) (eV)
%   Unk       : nb x nb x Nkx x Nky, eigenvectors (columns = |u_n>)
%   Enk       : Nkx x Nky x nb, eigenvalues (eV)
%   band_list : 参与求和/输出的带编号 (1 x Nb_sel)
%   dk_vecs   : 2 x 3, k-space step vectors (用于步长模长)
%   deltaE_reg: (optional) regularization threshold for |ΔE| (eV), default 1e-4
%
% Outputs:
%   Dk_core : Nkx x Nky x 2 x 2 x 2 x Nb_sel (single)
%   gk      : Nkx x Nky x 2 x 2 x Nb_sel     (single)
%
% Notes:
%   - 当前版本假设 Ndir=2 (x,y)
%   - gk 是在 band_list 子空间内对 m 求和的"截断"结果

    % ---- constants ----
    hbar_eVs = 6.582119569e-16;   % eV*s (用于 v = (1/ħ)∂H/∂k)

    [~, ~, Nkx, Nky] = size(Hk);
    Ndir   = size(dk_vecs, 1);
    Nb_sel = numel(band_list);

    if Ndir ~= 2
        error('Current version assumes Ndir=2 (2D x,y).');
    end

    if nargin < 6 || isempty(deltaE_reg)
        deltaE_reg = 1e-4;  % eV
    end

    % ---- dk step length for central difference ----
    dk_norm = sqrt(sum(dk_vecs.^2, 2));  % 2 x 1
    dkx = dk_norm(1);
    dky = dk_norm(2);

    % ---- output buffers ----
    Dk_core = zeros(Nkx, Nky, 2, 2, 2, Nb_sel, 'single');
    gk      = zeros(Nkx, Nky, 2, 2, Nb_sel,     'single');

    % ===== main loop (parfor on ix) =====
    parfor ix = 1:Nkx
        % local slices for this ix
        Dk_ix = zeros(Nky, 2, 2, 2, Nb_sel, 'single');
        gk_ix = zeros(Nky, 2, 2, Nb_sel,     'single');

        % periodic neighbors in x
        ixp = ix + 1; if ixp > Nkx, ixp = 1; end
        ixm = ix - 1; if ixm < 1,   ixm = Nkx; end

        for iy = 1:Nky
            % periodic neighbors in y
            iyp = iy + 1; if iyp > Nky, iyp = 1; end
            iym = iy - 1; if iym < 1,   iym = Nky; end

            % eigenvectors at (ix,iy)
            Uk_full = squeeze(Unk(:,:,ix,iy));    % nb x nb
            Uk      = Uk_full(:, band_list);      % nb x Nb_sel

            % energies for ΔE
            Ek_full = squeeze(Enk(ix,iy,:));      % nb x 1
            Ek      = Ek_full(band_list);         % Nb_sel x 1

            % ---- compute Hx, Hy by central difference (periodic) ----
            Hx = (Hk(:,:,ixp,iy) - Hk(:,:,ixm,iy)) / (2*dkx);
            Hy = (Hk(:,:,ix,iyp) - Hk(:,:,ix,iym)) / (2*dky);

            % ---- band-space M matrices ----
            Mx = Uk' * Hx * Uk;  % Nb_sel x Nb_sel
            My = Uk' * Hy * Uk;

            % ---- ΔE and denom^2 ----
            dE_nm = Ek - Ek.';                    % Nb_sel x Nb_sel
            small = abs(dE_nm) < deltaE_reg;
            dE_nm(small) = deltaE_reg .* sign(real(dE_nm(small)) + deltaE_reg);

            den2 = dE_nm.^2;
            den2(1:Nb_sel+1:end) = Inf;           % exclude m=n

            % ---- loop over n ----
            for in = 1:Nb_sel
                % velocities (real)
                vx_n = real(Mx(in,in)) / hbar_eVs;
                vy_n = real(My(in,in)) / hbar_eVs;

                % build metric components for this n
                Mx_row = Mx(in,:);
                My_row = My(in,:);
                Mx_col = Mx(:,in).';
                My_col = My(:,in).';
                den2_row = den2(in,:);

                gxx = real(sum( (Mx_row .* Mx_col) ./ den2_row ));
                gxy = real(sum( (Mx_row .* My_col) ./ den2_row ));
                gyx = real(sum( (My_row .* Mx_col) ./ den2_row ));
                gyy = real(sum( (My_row .* My_col) ./ den2_row ));

                % save g
                gk_ix(iy,1,1,in) = single(gxx);
                gk_ix(iy,1,2,in) = single(gxy);
                gk_ix(iy,2,1,in) = single(gyx);
                gk_ix(iy,2,2,in) = single(gyy);

                % D core
                v    = [vx_n, vy_n];
                gmat = [gxx,  gxy;
                        gyx,  gyy];

                for a = 1:2
                    for b = 1:2
                        for c = 1:2
                            D_core = v(a)*gmat(b,c) - v(b)*gmat(a,c);
                            Dk_ix(iy,a,b,c,in) = single(D_core);
                        end
                    end
                end
            end
        end

        % write back ix slice
        Dk_core(ix,:,:,:,:,:) = permute(Dk_ix, [1 2 3 4 5]); % Nky x 2 x 2 x 2 x Nb
        gk(ix,:,:,:,:)        = permute(gk_ix, [1 2 3 4]);   % Nky x 2 x 2 x Nb
    end
end


% function Dk_core = get_Dk_qmd_plainD_core( ...
%         Hk, Unk, Enk, band_list, dk_vecs, deltaE_reg)
% % GET_DK_QMD_PLAIND_CORE
% %   只计算并保存 k-resolved 的 plain-QMD 核:
% %
% %       D^{(n)}_{abc}(k) = v_a^n(k) g^{(n)}_{bc}(k) - v_b^n(k) g^{(n)}_{ac}(k)
% %
% %   其中
% %       g^{(n)}_{bc}(k) = Re Σ_{m≠n} M^b_{nm}(k) M^c_{mn}(k) / (ε_n-ε_m)^2
% %       M^b_{nm}(k) = <u_n(k)| ∂_{k_b} H(k) |u_m(k)>
% %       v_a^n(k)    = (1/ħ) Re M^a_{nn}(k)
% %
% % Inputs:
% %   Hk        : nb x nb x Nkx x Nky, H(k) (eV)
% %   Unk       : nb x nb x Nkx x Nky, eigenvectors (columns = |u_n>)
% %   Enk       : Nkx x Nky x nb, eigenvalues (eV)  (这里其实没直接用到，但保留接口方便以后扩展)
% %   band_list : 参与求和的带编号 (1 x Nb_sel)
% %   dk_vecs   : Ndir x 3, k-space step vectors (用于步长模长)
% %   weights   : (可选) Nkx x Nky, 用于后处理时积分；本函数不使用，但保留接口一致性
% %   deltaE_reg: (可选) 正则化能隙阈值 (eV), 默认 1e-4
% %
% % Output:
% %   Dk_core   : Nkx x Nky x 2 x 2 x 2 x Nb_sel (single)
% %              Dk_core(ix,iy,a,b,c,in) = D^{(n=in)}_{abc}(k)
% 
%     %#ok<NASGU>  % weights 当前不使用，仅为了接口兼容
% 
%     % ---- constants ----
%     hbar_eVs = 6.582119569e-16;   % eV*s  (用于 v = (1/ħ)∂H/∂k)
% 
%     [nb, ~, Nkx, Nky] = size(Hk);
%     Ndir   = size(dk_vecs, 1);
%     Nb_sel = numel(band_list);
% 
%     if Ndir ~= 2
%         error('Current version assumes Ndir=2 (2D x,y).');
%     end
% 
%     if nargin < 7 || isempty(deltaE_reg)
%         deltaE_reg = 1e-4;  % eV
%     end
% 
%     % ---- dk step length for central difference ----
%     dk_norm = sqrt(sum(dk_vecs.^2, 2));  % Ndir x 1
%     dkx = dk_norm(1);
%     dky = dk_norm(2);
% 
%     % ---- output buffer ----
%     Dk_core = zeros(Nkx, Nky, 2, 2, 2, Nb_sel, 'single');
% 
%     % ===== main loop (parfor on ix) =====
%     parfor ix = 1:Nkx
%         Dk_ix = zeros(Nky, 2, 2, 2, Nb_sel, 'single');
% 
%         % periodic neighbors in x
%         ixp = ix + 1; if ixp > Nkx, ixp = 1; end
%         ixm = ix - 1; if ixm < 1,   ixm = Nkx; end
% 
%         for iy = 1:Nky
%             % periodic neighbors in y
%             iyp = iy + 1; if iyp > Nky, iyp = 1; end
%             iym = iy - 1; if iym < 1,   iym = Nky; end
% 
%             % eigenvectors at (ix,iy)
%             Uk_full = squeeze(Unk(:,:,ix,iy));    % nb x nb
% 
%             % restrict to selected bands
%             Uk = Uk_full(:, band_list);           % nb x Nb_sel
% 
%             % ---- compute Hx, Hy by central difference (periodic) ----
%             Hx = (Hk(:,:,ixp,iy) - Hk(:,:,ixm,iy)) / (2*dkx);
%             Hy = (Hk(:,:,ix,iyp) - Hk(:,:,ix,iym)) / (2*dky);
% 
%             % ---- band-space M matrices: M^x = <u|Hx|u>, M^y = <u|Hy|u> ----
%             Mx = Uk' * Hx * Uk;  % Nb_sel x Nb_sel
%             My = Uk' * Hy * Uk;
% 
%             % ---- NOTE: dE_nm should use eigenvalues of selected bands ----
%             % We need energies at (ix,iy) to form (ε_n-ε_m)^2.
%             Ek_full = squeeze(Enk(ix,iy,:));      % nb x 1
%             Ek = Ek_full(band_list);              % Nb_sel x 1
% 
%             dE_nm = Ek - Ek.';                    % Nb_sel x Nb_sel
% 
%             % regularize small gaps (including diagonal)
%             small = abs(dE_nm) < deltaE_reg;
%             dE_nm(small) = deltaE_reg .* sign(real(dE_nm(small)) + deltaE_reg);
% 
%             den2 = dE_nm.^2;
%             den2(1:Nb_sel+1:end) = Inf;           % exclude m=n for metric sum
% 
%             % ---- loop over bands n ----
%             for in = 1:Nb_sel
%                 % velocities v_x^n, v_y^n (real scalars)
%                 vx_n = real(Mx(in,in)) / hbar_eVs;
%                 vy_n = real(My(in,in)) / hbar_eVs;
% 
%                 % rows/cols for g
%                 Mx_row = Mx(in,:);
%                 My_row = My(in,:);
%                 Mx_col = Mx(:,in).';
%                 My_col = My(:,in).';
% 
%                 den2_row = den2(in,:);
% 
%                 gxx = real(sum( (Mx_row .* Mx_col) ./ den2_row ));
%                 gxy = real(sum( (Mx_row .* My_col) ./ den2_row ));
%                 gyx = real(sum( (My_row .* Mx_col) ./ den2_row ));
%                 gyy = real(sum( (My_row .* My_col) ./ den2_row ));
% 
%                 v    = [vx_n, vy_n];
%                 gmat = [gxx,  gxy;
%                         gyx,  gyy];
% 
%                 for a = 1:2
%                     for b = 1:2
%                         for c = 1:2
%                             D_core = v(a)*gmat(b,c) - v(b)*gmat(a,c);
%                             Dk_ix(iy,a,b,c,in) = single(D_core);
%                         end
%                     end
%                 end
%             end
%         end
% 
%         % write back this ix slice
%         Dk_core(ix,:,:,:,:,:) = permute(Dk_ix, [1 2 3 4 5]); % Nky x 2 x 2 x 2 x Nb
%     end
% end
