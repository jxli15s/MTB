function [Dk_abc_n, D_abc_n] = ...
    get_D_metric_kresolved_from_gv( ...
        gk, vk, Enk, band_list, Ef_list, AreaBZ, eta, weights)
% D_METRIC_KRESOLVED_FROM_GV
%   利用已经计算好的 vk 和 gk，在 2D k-mesh 上给出
%   每条能带 n 的 k-resolved metric dipole integrand
%
%     D_{a;bc}^{(n)}(k;E_F)
%       = (v_a g_bc - v_b g_ac) δ(ε_n(k) - E_F)
%
%   以及对 k 积分后的
%
%     D_{a;bc}^{(n)}(E_F)
%       = sum_k D_{a;bc}^{(n)}(k;E_F) * ΔS / (2π)^2.
%
% 输入：
%   gk      : Nkx x Nky x Ndir x Ndir x Nb_sel
%             gk(ix,iy,a,b,ib) = g_{ab}^{(band_list(ib))}(k)
%   vk      : Nkx x Nky x Ndir       x Nb_sel
%             vk(ix,iy,a,ib)   = v_a^{(band_list(ib))}(k)
%   Enk     : Nkx x Nky x nb     —— 所有能带能量
%   band_list : 参与计算的带编号，长度 Nb_sel
%   Ef_list   : 1 x N_EF 的 Fermi 能量列表
%   AreaBZ    : 2D Brillouin 区面积
%   eta       : δ(ε−E_F) 的洛伦兹展宽参数
%   weights   : Nkx x Nky 的 k 点权重(面积元 ΔS)；
%               若为空[]则自动设为 AreaBZ/Nk。
%
% 输出：
%   Dk_abc_n : Nkx x Nky x Ndir x Ndir x Ndir x Nb_sel x N_EF
%              Dk_abc_n(ix,iy,a,b,c,ib,ief)
%              = D_{a;bc}^{(n)}(k;E_F_ief) 的 integrand
%
%   D_abc_n  : Ndir x Ndir x Ndir x Nb_sel x N_EF
%              D_abc_n(a,b,c,ib,ief)
%              = 对 k 积分后的 D_{a;bc}^{(n)}(E_F_ief)

    [Nkx, Nky, Ndir, ~, Nb_sel] = size(gk);
    N_EF = numel(Ef_list);

    % 只取参与的带的能量：Nkx x Nky x Nb_sel
    Enk_sel = Enk(:,:,band_list);

    % k 权重，如果没给就均匀
    if nargin < 8 || isempty(weights)
        Nk_tot  = Nkx * Nky;
        weights = (AreaBZ / Nk_tot) * ones(Nkx, Nky);
    end

    % 结果：k-resolved & k-integrated
    Dk_abc_n = zeros(Nkx, Nky, Ndir, Ndir, Ndir, Nb_sel, N_EF);
    D_abc_n  = zeros(Ndir, Ndir, Ndir, Nb_sel, N_EF);

    % ========= 主循环 =========
    for ief = 1:N_EF
        Ef = Ef_list(ief);

        % δ(ε_n - E_F)，形状：Nkx x Nky x Nb_sel
        dE     = Enk_sel - Ef;
        deltaF = (1/pi) * eta ./ (dE.^2 + eta^2);

        % 遍历每一条带 n
        for ib = 1:Nb_sel
            deltaF_band = deltaF(:,:,ib);     % Nkx x Nky

            for a = 1:Ndir           % 电流方向
                v_a = vk(:,:,a,ib);

                for b = 1:Ndir       % 第一个场方向
                    v_b = vk(:,:,b,ib);

                    for c = 1:Ndir   % 第二个场方向

                        g_bc = gk(:,:,b,c,ib);
                        g_ac = gk(:,:,a,c,ib);

                        % k-resolved integrand:
                        integrand = (v_a .* g_bc - v_b .* g_ac) .* deltaF_band;
                        Dk_abc_n(:,:,a,b,c,ib,ief) = integrand;

                        % 对 k 积分：∑_k integrand * ΔS / (2π)^2
                        D_val = sum(integrand .* weights, 'all') / (2*pi)^2;
                        D_abc_n(a,b,c,ib,ief) = D_abc_n(a,b,c,ib,ief) + D_val;
                    end
                end
            end
        end
    end
end
