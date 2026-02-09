function sigma_abc = get_sigma_quantum_metric_dipole( ...
        Hk, Unk, Enk, band_list, dk_vecs, Ef_list, ...
        AreaBZ, eta, weights, deltaE_reg, T_K)
% SIGMA_ABC_FROM_HMESH_PERIODIC
%   使用 mesh 上预先给出的 H(k)，在 2D BZ 上采用周期性边界的中心差分
%   H_i = ∂H/∂k_i，计算
%
%   σ_{abc}(E_F) =
%     (2 e^3 / ħ^2) ∑_{n,m≠n,k}
%       Re{ [ v_a^n M^b_{nm} M^c_{mn}
%             - v_b^n M^a_{nm} M^c_{mn} ]
%           / (ε_n - ε_m)^3 }
%       δ(ε_n - E_F) d^2k / (2π)^2
%
% 输入:
%   Hk        : nb x nb x Nkx x Nky, H(k) (单位 eV).
%   Unk       : nb x nb x Nkx x Nky, 本征矢矩阵 (列为 |u_n>).
%   Enk       : Nkx x Nky x nb, 本征值 (eV).
%   band_list : 参与求和的带编号, 长度 Nb_sel.
%   dk_vecs   : Ndir x 3, 每个方向在 k 空间的步长向量,
%               e.g. [b1/knum; b2/knum] for 2D.
%   Ef_list   : 1 x N_EF 的 Fermi 能列表 (eV).
%   AreaBZ    : Brillouin 区总面积 (与 k 单位匹配, 如 1/Å^2).
%   eta       : δ(ε−E_F) 展宽参数 (eV).
%   weights   : Nkx x Nky 的 k 权重; 若 [] 则默认均匀 AreaBZ/(Nkx*Nky).
%   deltaE_reg: 最小能隙正则化阈值 (eV), 如 1e-4.
%
% 输出:
%   sigma_abc : Ndir x Ndir x Ndir x N_EF
%               sigma_abc(a,b,c,ief) = 对应 E_F 的 σ_{abc}(E_F).

    % ==== 常数 (用 eV·s 的 ħ, 方便和 H 一致) ====
    e_charge = 1.602176634e-19;   % C
    hbar_eVs = 6.582119569e-16;   % eV·s

    % pref = 2* e_charge^0
    pref = 2*e_charge^1;

    [nb, ~, Nkx, Nky] = size(Hk);
    Ndir   = size(dk_vecs, 1);
    Nb_sel = numel(band_list);
    N_EF   = numel(Ef_list);

    dk_norm = sqrt(sum(dk_vecs.^2, 2));   % 每个方向步长的模长

    % 只取参与能带的能量: Nkx x Nky x Nb_sel
    Enk_sel = Enk(:,:, band_list);

    % ------------ k 权重 ------------
    if nargin < 9 || isempty(weights)
        Nk_tot  = Nkx * Nky;
        weights = (AreaBZ / Nk_tot) * ones(Nkx, Nky);
    end

    if nargin < 10 || isempty(deltaE_reg)
        deltaE_reg = 1e-4;   % eV
    end

    % σ 的几何部分(先不乘 prefactor): sum_{n,m,k} ... δ(...)
    geom_sigma = zeros(Ndir, Ndir, Ndir, N_EF);

    % ========== 主循环: 扫描所有 k 点 (含边界, 周期性差分) ==========
    parfor ix = 1:Nkx
        fprintf("Processing on ix= %d of total %d Nkx\n",ix,Nkx)
        tic;
        % 每个 ix 用一个局部累加器
        geom_sigma_loc = zeros(Ndir, Ndir, Ndir, N_EF);

        % 周期边界下的 x 邻居索引
        ixp = ix + 1; if ixp > Nkx, ixp = 1; end
        ixm = ix - 1; if ixm < 1,   ixm = Nkx; end

        for iy = 1:Nky
            % 周期边界下的 y 邻居索引
            iyp = iy + 1; if iyp > Nky, iyp = 1; end
            iym = iy - 1; if iym < 1,   iym = Nky; end

            % 本征值 / 本征矢 (只取 band_list 子空间)
            Uk_full = squeeze(Unk(:,:,ix,iy));    % nb x nb
            Ek_full = squeeze(Enk(ix,iy,:));      % nb x 1
            Uk = Uk_full(:, band_list);           % nb x Nb_sel
            Ek = Ek_full(band_list);              % Nb_sel x 1

            % ----- 计算 H_i(=∂H/∂k_i) 和 M^i_nm = <n|H_i|m> -----
            % M_all(n,m,i)
            M_all = zeros(Nb_sel, Nb_sel, Ndir);

            for id = 1:Ndir
                if id == 1
                    % 沿 "x" 方向的中心差分: (H_{ix+1,iy} - H_{ix-1,iy})/(2|dk_x|)
                    H_plus  = Hk(:,:,ixp,iy);
                    H_minus = Hk(:,:,ixm,iy);
                elseif id == 2
                    % 沿 "y" 方向的中心差分
                    H_plus  = Hk(:,:,ix,iyp);
                    H_minus = Hk(:,:,ix,iym);
                else
                    error('当前版本只实现 Ndir=2 (2D)。');
                end

                H_i = (H_plus - H_minus) ./ (2 * dk_norm(id));  % 近似 ∂H/∂k_i

                % 转到带空间: M^i_nm = <u_n|H_i|u_m>
                M_all(:,:,id) = Uk' * H_i * Uk;   % Nb_sel x Nb_sel
            end

            % ----- 能量差矩阵 ΔE_nm = ε_n - ε_m -----
            Ec = Ek;          % (Nb_sel x 1)
            Er = Ek.';        % (1 x Nb_sel)
            dE_nm = Ec - Er;  % ΔE_nm(n,m) = ε_n - ε_m

            % 正则化避免 0 (包括 n=m，会被正则化，但后面只用 m≠n)
            small = abs(dE_nm) < deltaE_reg;
            dE_nm(small) = deltaE_reg .* sign(real(dE_nm(small)) + deltaE_reg);

            % ----- 对每个带 n, 再对 m≠n 求和 -----
            for in = 1:Nb_sel
                En = Ek(in);          % 当前带能量

                % % δ(ε_n - E_F) -> Lorentz 展宽: 1 x N_EF
                % dE_EF   = En - Ef_list;                         % 1 x N_EF
                % delta_n = (1/pi) * eta ./ (dE_EF.^2 + eta^2);   % 1 x N_EF

                kB_eV = 8.617333262e-5;   % eV/K
                if T_K == 0
                    dE_EF   = En - Ef_list;                       % 1 x N_EF
                    delta_n = (1/pi) * eta ./ (dE_EF.^2 + eta^2); % 1 x N_EF
                else
                    kBT = kB_eV * T_K;
                    % -df/dE = 1/(4 kBT) * sech^2(x)
                    x   = (En - Ef_list) ./ (2*kBT);               % 1 x N_EF
                    delta_n = (1./(4*kBT)) ./ cosh(x).^2;          % 1 x N_EF
                end

                % 带内速度 v_a^n = (1/ħ) Re M^a_nn
                v_n = zeros(1, Ndir);
                for a = 1:Ndir
                    M_a = M_all(:,:,a);
                    v_n(a) = real(M_a(in,in))/hbar_eVs;      % 速度应为实数
                end

                % 对所有 m≠n 求和
                for im = 1:Nb_sel
                    if im == in, continue; end

                    dE    = dE_nm(in, im);     % ε_n - ε_m
                    dE3   = dE^3;

                    for a = 1:Ndir
                        M_a = M_all(:,:,a);

                        for b = 1:Ndir
                            M_b = M_all(:,:,b);

                            for c = 1:Ndir
                                M_c = M_all(:,:,c);

                                % M^b_{nm}, M^c_{mn}, M^a_{nm}
                                M_b_nm = M_b(in,im);
                                M_c_mn = M_c(im,in);
                                M_a_nm = M_a(in,im);

                                % 核心 integrand：
                                % [ v_a M^b_nm M^c_mn - v_b M^a_nm M^c_mn ] / (ΔE^3)
                                core = ( v_n(a) * M_b_nm * M_c_mn ...
                                    - v_n(b) * M_a_nm * M_c_mn ) / dE3;

                                % 乘上 δ(ε_n - E_F) 并取实部 -> 1 x N_EF
                                contrib_vec = real(core) .* delta_n;  % 1 x N_EF

                                % 加到几何部分（这里已经乘了该点的 k 权重）
                                w = weights(ix,iy) / (2*pi)^2;
                                % 这里累加到局部 geom_sigma_loc
                                geom_sigma_loc(a,b,c,:) = geom_sigma_loc(a,b,c,:) ...
                                    + reshape(contrib_vec * w, [1,1,1,N_EF]);
                            end
                        end
                    end
                end
            end
        end
        % parfor reduction: geom_sigma = geom_sigma + geom_sigma_loc
        geom_sigma = geom_sigma + geom_sigma_loc;
        toc;
    end

    % 乘上整体 prefactor 得到 σ_{abc}(E_F)
    sigma_abc = pref * geom_sigma;
end
