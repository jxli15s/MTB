function [gk, Qk, Fk, vk] = quantum_geometry_general_plane( ...
        obj, Kx, Ky, Kz, Unk, Enk, band_list, dk_list, delta)
% PLANE_QUANTUM_GEOMETRY
%   在二维 k-mesh 上，对若干条能带计算量子度规 g_ij、
%   量子几何张量 Q_ij 和 Berry 曲率张量 F_ij。
%
% 输入：
%   obj       : 有 get_hk(k) 方法的 TB 对象（和 quantum_geometry_general 里一致）
%   Kx,Ky,Kz  : Nkx x Nky 的 k 网格
%   Unk       : nb x nb x Nkx x Nky, 每个 k 点的本征矢矩阵
%   Enk       : Nkx x Nky x nb, 每个 k 点的本征值
%   band_list : 需要计算的能带编号数组，例如 [1 2 3]
%   dk_list   : Ndir x D 的小位移方向（和 quantum_geometry_general 一样）
%   delta     : 规避小能隙的阈值
%
% 输出：
%   gk  : Nkx x Nky x Ndir x Ndir x Nb_sel
%         其中 Nb_sel = numel(band_list)，gk(ix,iy,:,:,ib) 就是该点该带的 g_ij
%   Qk  : 同上，但存量子几何张量 Q_ij
%   Fk  : 同上，但存 Berry 曲率张量 F_ij = 2 Im Q_ij
%
% 使用示例（2D 系统）：
%   b1 = obj.b(1,:);  b2 = obj.b(2,:);
%   dk_list = [b1/knum; b2/knum];   % 对应 kx, ky 方向
%   [gk,Qk,Fk] = plane_quantum_geometry(g,Kx,Ky,Kz,Unk,Enk,[1,2],dk_list,1e-6);
%   gxx_band1 = squeeze(gk(:,:,1,1,1));
%   gyy_band1 = squeeze(gk(:,:,2,2,1));
%   gxy_band1 = squeeze(gk(:,:,1,2,1));
%   Omega_z_band1 = squeeze(Fk(:,:,1,2,1));   % F_{xy}
%

    [Nkx, Nky] = size(Kx);
    Ndir       = size(dk_list, 1);
    nb_sel     = numel(band_list);
    quantum_geometry_general=@MTB.ham.quantum_geometry_general_k;
    % 预分配：Nkx x Nky x Ndir x Ndir x Nb_sel
    gk = zeros(Nkx, Nky, Ndir, Ndir, nb_sel);
    Qk = zeros(Nkx, Nky, Ndir, Ndir, nb_sel);
    Fk = zeros(Nkx, Nky, Ndir, Ndir, nb_sel);
    vk = zeros(Nkx, Nky, Ndir,       nb_sel);  % 每带每方向的 v_i

    % 主循环：遍历 k 网格
    parfor ix = 1:Nkx
        for iy = 1:Nky

            % 当前 k 点
            k = [Kx(ix,iy), Ky(ix,iy), Kz(ix,iy)];

            % 当前 k 点的本征矢、本征值
            % 假设 Enk(ix,iy,:) 是 nb x 1 的能量，Unk(:,:,ix,iy) 是 nb x nb 的本征矢矩阵
            Ek = squeeze(Enk(ix,iy,:));          % nb x 1
            Uk = squeeze(Unk(:,:,ix,iy));        % nb x nb

            % 对 band_list 中的每一条带都算一遍量子几何
            for ib = 1:nb_sel
                n = band_list(ib);

                [g_loc, Q_loc, F_loc, v_dirs] = quantum_geometry_general( ...
                    obj, k, dk_list, Uk, Ek, n, delta);

                gk(ix,iy,:,:,ib) = g_loc;
                Qk(ix,iy,:,:,ib) = Q_loc;
                Fk(ix,iy,:,:,ib) = F_loc;
                vk(ix,iy,:,ib)   = v_dirs;
            end
        end
    end
end
