classdef geometry < handle
    properties
        name;
        a;
        a2;
        b;
        b2;
        atoms;
        hopr;
        hopr2;
        ham;
        wpos;
        iniham;
        sublattice;
        orbnum_list;
        suborbidx;
        Rcart;
    end

    methods
        function obj = geometry(name)
            obj.name = name;
            obj.a = [];
            obj.a2 = [];
            obj.b = [];
            obj.b2 = [];
            obj.atoms = [];
            obj.wpos = [];
            obj.hopr = [];
            obj.hopr2 = [];
            obj.ham = [];
            obj.iniham = [];
            obj.sublattice=[];
            obj.suborbidx=[];
            obj.orbnum_list=[];
            obj.Rcart=[];
        end

        function Umatrix=MillerIndicestoumatrix(obj,MillerIndices)
            Umatrix=zeros(3,3);
            rmiller=gcd(sym(MillerIndices));
            h=MillerIndices(1)/rmiller;
            k=MillerIndices(2)/rmiller;
            l=MillerIndices(3)/rmiller;
            vector_on_hkl_surface=[];
            for i1 = -5:5
                for i2 = -5:5
                    for i3 = -5:5
                        if i1==0 && i2==0 && i3==0
                            continue
                        elseif abs(dot([i1,i2,i3],MillerIndices)) < eps
                            vector_on_hkl_surface=[vector_on_hkl_surface;[i1,i2,i3]];
                        end
                    end
                end
            end
            length(vector_on_hkl_surface);
            smallest_area = 99999999.0;
            for i1 = 1 : length(vector_on_hkl_surface)-1
                for i2 = i1+1 : length(vector_on_hkl_surface)
                    R1 = vector_on_hkl_surface(i1,:)*obj.a;
                    R2 = vector_on_hkl_surface(i2,:)*obj.a;
                    R3 = cross(R1,R2);
                    area = norm(R3);
                    if area < 0.01
                        continue
                    elseif area < smallest_area
                        smallest_area = area;
                    end
                end
            end

            largestangle=0.0;
            for i1 = 1 : length(vector_on_hkl_surface)-1
                for i2 = i1+1 : length(vector_on_hkl_surface)
                    R1 = vector_on_hkl_surface(i1,:)*obj.a;
                    R2 = vector_on_hkl_surface(i2,:)*obj.a;
                    R3 = cross(R1,R2);
                    area = norm(R3);
                    % angle = acos(dot(R1,R2)/norm(R1)/norm(R2))
                    angle = acos(dot(R1,R2)/norm(R1)/norm(R2));
                    if angle > pi/2.0
                        angle = abs(angle-pi);
                    end

                    if area - smallest_area < 10^-6
                        if angle > largestangle
                            largestangle = angle;
                        end
                    end
                end
            end

            flag = false;
            for i1 = 1 : length(vector_on_hkl_surface)-1
                for i2 = i1+1 : length(vector_on_hkl_surface)
                    R1 = vector_on_hkl_surface(i1,:)*obj.a;
                    R2 = vector_on_hkl_surface(i2,:)*obj.a;
                    R3 = cross(R1,R2);
                    area = norm(R3);
                    angle = acos(dot(R1,R2)/norm(R1)/norm(R2));
                    if angle > pi/2.0
                        angle = abs(angle-pi);
                    end


                    if abs(area - smallest_area) < 10^-6 && abs(angle - largestangle) < 10^-6
                        Umatrix(1,:) = vector_on_hkl_surface(i1,:);
                        Umatrix(2,:) = vector_on_hkl_surface(i2,:);
                        flag = true;
                        break;
                    end
                end
                if flag
                    break;
                end
            end

            smallest_volume = 999999999.0;
            R1 = Umatrix(1,:)*obj.a;
            R2 = Umatrix(2,:)*obj.a;
            for i1 = -5:5
                for i2 = -5:5
                    for i3 = -5:5
                        if i1==0 && i2==0 && i3==0
                            continue
                        end
                        R3 = [i1,i2,i3]*obj.a;
                        cell_volume = abs(dot(R3,cross(R1,R2)));
                        if cell_volume < 10^-6
                            continue
                        end
                        if cell_volume < smallest_volume
                            smallest_volume = cell_volume;
                        end
                    end
                end
            end

            smallest_length = 99999999.0;
            for i1 = -5:5
                for i2 = -5:5
                    for i3 = -5:5
                        if i1==0 && i2==0 && i3==0
                            continue
                        else
                            R3 = [i1,i2,i3]*obj.a;
                            cell_volume = abs(dot(R3,cross(R1,R2)));
                            if abs(cell_volume - smallest_volume) < 10^-6 && ...
                                    norm(R3) < smallest_length
                                smallest_length = norm(R3);
                            end
                        end
                    end
                end
            end

            flag = false;
            for i1 = -5:5
                for i2 = -5:5
                    for i3 = -5:5
                        if i1==0 && i2==0 && i3==0
                            continue
                        else
                            R3 = [i1,i2,i3]*obj.a;
                            cell_volume = abs(dot(R3,cross(R1,R2)));
                            if abs(cell_volume - smallest_volume) < 10^-6 && ...
                                    abs(norm(R3) - smallest_length) < 10^-6
                                Umatrix(3,:) = [i1,i2,i3];
                                flag = true;
                                break
                            end
                        end
                    end
                    if flag
                        break;
                    end
                end
                if flag
                    break;
                end
            end


            cell_volume_new=det(Umatrix*obj.a);
            if (cell_volume_new<0)
                Umatrix(3,:)=-Umatrix(3,:);
            end
            a_new=Umatrix*obj.a;
            if abs(det(a_new)-det(obj.a)) < 10^-6
                obj.a=a_new;
                obj.b=inv(obj.a')*2*pi;
                obj.atoms=obj.atoms*inv(Umatrix);
                obj.hopr2 = obj.hopr*inv(Umatrix);
            end
            % Umatrix=zeros(3,3);
            % Umatrix=largestangle
        end

        function Urot = surfab(obj)
            Urot=eye(3);
            Urot(1,:) = obj.a(1,:)/norm(obj.a(1,:));
            Urot(3,:) = cross(obj.a(1,:),obj.a(2,:))/norm(cross(obj.a(1,:),obj.a(2,:)));
            Urot(2,:) = cross(Urot(3,:),Urot(1,:))/norm(cross(Urot(3,:),Urot(1,:)));
            obj.a2 = obj.a*inv(Urot);
            obj.a2 = obj.a2(1:2,1:2);
            obj.b2 = inv(obj.a2')*2*pi;
        end

        function readwpos(obj,filename)
            obj.wpos=table2array(readtable(filename))
        end

        % function hk=get_hk(obj,kpoint)
        %     ham1=obj.ham;
        %     hopring=obj.hopr;
        %     hk=zeros(size(ham1,1),size(ham1,1));
        %     nrpts=size(obj.ham,3);
        %     ar=obj.a;
        %     for j=1:nrpts
        %         temp=ham1(:,:,j)*...
        %             exp(1j*dot(kpoint,(hopring(j,:)*ar)));
        %         hk=hk+temp;
        %     end
        % end
        function [E,Psik]=solve_kpoint(obj,k)
            %slove_kpoint: this function give the energy level and eigenstate at the given
            %k-point.
            HK=obj.get_hk(k);
            [V,D]=eig(HK);
            [E,ind]=sort(diag(real(D)));
            % [E,ind]=sort(diag(D));
            E=real(E);
            % E=real(E);
            Psik=V(:,ind);
        end

        function hk=get_hk(obj,kpoint)
            hamiltonian=reshape(obj.ham,[],size(obj.hopr,1));
            hk=hamiltonian*exp(1j*((obj.hopr*obj.a)*kpoint'));
            hk=reshape(hk,size(obj.ham,1),[]);
            hk=(hk+hk')/2.0;
        end


        % function hk=get_hk_bac(obj,kpoint)
        %     hk=zeros(size(obj.ham,1),size(obj.ham,1));
        %     nrpts=size(obj.ham,3);
        %     for j=1:nrpts
        %         temp=obj.ham(:,:,j)*...
        %             exp(1j*dot(kpoint,(obj.hopr(j,:)*obj.a)));
        %         hk=hk+temp;
        %     end
        % end
        %
        %  function hk=get_hk(obj,kpoint)
        %     nbands=size(obj.ham,1);
        %     nrpts=size(obj.ham,3);
        %     phase=exp(1j*((obj.hopr*obj.a)*kpoint'));
        %     phase=kron(phase,ones(nbands^2,1));
        %     phase=reshape(phase,nbands,nbands,nrpts);
        %     hk=obj.ham.*phase;
        %     hk=sum(hk,3);
        % end



        function [Kx,Ky] = get_Slab2Dkmesh(obj,kxline,kyline,knum)
            x_frac=linspace(kxline(1),kxline(2),knum+1);
            y_frac=linspace(kyline(1),kyline(2),knum+1);
            x_frac=x_frac(1:knum);
            y_frac=y_frac(1:knum);
            [X_frac,Y_frac]=meshgrid(x_frac,y_frac);
            Kx=X_frac.*obj.b2(1,1)+Y_frac.*obj.b2(2,1);
            Ky=X_frac.*obj.b2(1,2)+Y_frac.*obj.b2(2,2);
        end

        function [Kx,Ky,Kz] = get_Bulk2Dkmesh(obj,kxline,kyline,knum)
            x_frac=linspace(kxline(1),kxline(2),knum+1);
            y_frac=linspace(kyline(1),kyline(2),knum+1);
            x_frac=x_frac(1:knum);
            y_frac=y_frac(1:knum);
            [X_frac,Y_frac]=meshgrid(x_frac,y_frac);
            Kx=X_frac.*obj.b(1,1)+Y_frac.*obj.b(2,1);
            Ky=X_frac.*obj.b(1,2)+Y_frac.*obj.b(2,2);
            Kz=X_frac.*obj.b(1,3)+Y_frac.*obj.b(2,3);
        end

        function obj=onsite_modify(obj,E)
            index=find(ismember(obj.hopr,[0,0,0],'rows'));
            obj.ham=obj.iniham+0;
            Ef=obj.ham(:,:,index)+E;
            obj.ham(:,:,index)=Ef(:,:);
        end

        function obj=add_zeeman(obj,E)
            index=find(ismember(obj.hopr,[0,0,0],'rows'));
            Ef=obj.ham(:,:,index)+E;
            obj.ham(:,:,index)=Ef(:,:);
        end


        function obj=offsite_modify(obj,hopr,ham)
            index=find(ismember(obj.hopr,hopr,'rows'));
            ham_new=obj.ham(:,:,index)+ham;
            obj.ham(:,:,index)=ham_new(:,:);
        end


        function get_suborbidx(obj)
            obj.suborbidx=zeros(size(obj.atoms,1),2);
            end_idx=0;
            for i=1:size(obj.atoms,1)
                start_idx=end_idx+1;
                end_idx=start_idx+obj.orbnum_list(i)-1;
                obj.suborbidx(i,:)=[start_idx,end_idx];
            end
        end

        function lattice_plot(obj)
            lattice=obj.a;
            sublattice=obj.atoms;
            figure('Color','white')
            drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} ) ;
            drawArrow([0,lattice(1,1)],[0,lattice(1,2)],'linewidth',2,'color','k')
            hold on
            drawArrow([0,lattice(2,1)],[0,lattice(2,2)],'linewidth',2,'color','k')
            hold on
            % sublattice(:,1:2)=sublattice(:,1:2)*(lattice);
            sublattice=sublattice*lattice;


            for i=1:size(sublattice,1)
                if sublattice(i,end)==1

                    plot(sublattice(i,1),sublattice(i,2),'.','Color','r','MarkerSize',40,'MarkerFaceColor','k')
                elseif sublattice(i,end)==2
                    plot(sublattice(i,1),sublattice(i,2),'.','Color','b','MarkerSize',40,'MarkerFaceColor','k')
                else
                    plot(sublattice(i,1),sublattice(i,2),'.','Color','#038103','MarkerSize',40,'MarkerFaceColor','k')
                end
                text(sublattice(i,1)+0.05,sublattice(i,2)+0.05,num2str(i))
            end
            axis equal
            axis off
        end

        %%
        function [bz_vertices, bz_ridges, bz_facets] = get_brillouin_zone_3d(obj)
            %GET_BRILLOUIN_ZONE_3D_LIKE_PYTHON
            % MATLAB 版：用 Voronoi 得到第一布里渊区（WS 单胞），并返回
            %   bz_vertices : (Nv x 3) 该区域用到的所有有限顶点坐标
            %   bz_ridges   : 1xNfaces cell，每个元素是该面的边界折线 Mx3（闭合但不重复首尾）
            %   bz_facets   : 1xNfaces cell，与 ridges 同序的面顶点集合（未排序）
            %
            % 说明：
            % - 与 Python 代码逻辑等价：只取"原点点"的 Voronoi 区域。
            % - MATLAB 的 voronoin 不直接给 ridge_points；我们依据"原点-邻点"的
            %   中垂面，从该区域的顶点里筛出面并排序成折线。
            %
            % 输入:
            %   cellB : 3x3 倒易基矢；列为 b1,b2,b3（若发现像行输入，会自动转置）

            % --- 输入与列向量检查 ---
            B = obj.b;
            assert(all(size(B) == [3 3]), 'cell must be 3x3');
            % 若用户按行给矢量，自动转为列
            if sum(vecnorm(B,2,1)) < sum(vecnorm(B,2,2).')
                B = B.';
            end

            % --- 生成 [-1,0,1]^3 倒易点阵（与原 Python 一致） ---
            [I,J,K] = ndgrid(-1:1, -1:1, -1:1);
            multi  = [I(:), J(:), K(:)];       % 27x3
            points = multi * B.';              % 27x3

            % 原点索引：离 (0,0,0) 最近的那个（等价于 Python 里固定的 13/14th）
            [~, origin_idx] = min(vecnorm(points,2,2));

            % --- Voronoi 分解 ---
            % V: 全部顶点坐标；C{i}: 第 i 个点的 Voronoi 区域顶点索引（1 表示无穷远）
            [V, C] = voronoin(points);

            % 原点对应区域（第一 BZ），去掉无穷远顶点
            reg_idx = C{origin_idx};
            reg_idx = reg_idx(reg_idx ~= 1);
            assert(~isempty(reg_idx), 'Region at origin has no finite vertices.');
            bz_vertices = V(reg_idx, :);       % Nv x 3
            Vreg = bz_vertices;

            % --- 构造与原点邻点对应的面（facets）与其边界折线（ridges） ---
            M = size(points,1);
            tol_plane = 1e-9;    % 判断点落在平面上的容差
            tol_sort  = 1e-12;   % 面内排序时的符号容差

            bz_facets = {};
            bz_ridges = {};
            used_keys = containers.Map('KeyType','char','ValueType','logical'); % 去重

            for pid = 1:M
                if pid == origin_idx, continue; end
                p = points(pid, :);                 % 邻点矢量 G
                n = p;                              % 中垂面法向量
                d = 0.5 * dot(p, p);                % 平面方程 n·x = d

                vals = Vreg * n.' - d;              % 各顶点到平面的代数距离
                on_plane = abs(vals) < tol_plane;
                idx_on   = find(on_plane);
                if numel(idx_on) < 3
                    continue;                       % 不是一个有效面
                end

                % --- 用规范化 (n,d) 作为键，避免正负法向重复 ---
                n_unit = n / max(norm(n), eps);
                d_unit = d / max(norm(n), eps);
                if n_unit(1) < -tol_sort || (abs(n_unit(1))<=tol_sort && ...
                        (n_unit(2) < -tol_sort || (abs(n_unit(2))<=tol_sort && n_unit(3) < -tol_sort)))
                    n_unit = -n_unit; d_unit = -d_unit;
                end
                key = sprintf('%.12f_%.12f_%.12f_%.12f', n_unit(1), n_unit(2), n_unit(3), d_unit);
                if isKey(used_keys, key)
                    continue;                       % 同一几何平面已记录
                end
                used_keys(key) = true;

                % --- 该面的顶点集合（未排序） ---
                facet_pts = Vreg(idx_on, :);
                bz_facets{end+1} = facet_pts; %#ok<AGROW>

                % --- 得到该面边界的"闭合折线"：把点投影到面内 2D，然后按极角排序 ---
                % 构造面内正交基（用 null(n) 即可）
                T = null(n_unit);                   % 3x2
                if size(T,2) ~= 2
                    % 退化情形，用 SVD 兜底
                    [U,~,~] = svd(eye(3) - n_unit.'*n_unit);
                    T = U(:,1:2);
                end
                c0 = mean(facet_pts, 1);            % 面内参考中心
                X2 = (facet_pts - c0) * T;          % 投影到 2D
                ang = atan2(X2(:,2), X2(:,1));
                [~, ord] = sort(ang);
                ridge_loop = facet_pts(ord, :);     % 闭合（首尾不重复）

                bz_ridges{end+1} = ridge_loop;      % 存入折线
            end
        end

    end
end



