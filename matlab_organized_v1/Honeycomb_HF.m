clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the geomtery information                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = MTB.geometry("honeycomb");
a1=[3/2, sqrt(3)/2, 0.0];
a2=[-3/2, sqrt(3)/2, 0.0];
a3=[0.0, 0.0, 10];
omega=dot(a1,cross(a2,a3));
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
b=[b1;b2;b3];
g.a = [a1;a2;a3];
g.b = b;
g.atoms=[-1/6,  1/6, 0;...
          1/6, -1/6, 0];
g.wpos=g.atoms*g.a;
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the real space Ham                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_matrix_1 = find_neighbor_data(g.wpos(1:end,:),g.a, 4,2);
%%
hopr = unique(result_matrix_1(:,3:5), 'rows');
nbands = 4;
ham = zeros(nbands,nbands,size(hopr,1));
t1=1;
for i =1:size(result_matrix_1,1)
    [~,raw_index]=ismember(result_matrix_1(i,3:5),hopr,"rows");
    orbital_1 = result_matrix_1(i,1);
    orbital_2 = result_matrix_1(i,2);
    ham(orbital_1*2-1:orbital_1*2,orbital_2*2-1:orbital_2*2,raw_index) = [1,0;0,1]*t1;
end
g.hopr = hopr;
g.ham = ham;
g.wpos=kron(g.wpos,ones(2,1));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [1/3,1/3,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
g.iniham=g.ham+0;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          Run Hartree Fock for AFM with Different U                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%how the  Band structure evolves in an interacting Graphene from a AFM guess
Us=0:0.5:4;
mz=[];
for i=1:size(Us,2)
g.ham=g.iniham; 
stepmax=1000;
minstepmax=11;
step=1;
knum=101;
Electric_field_in_evpA=0;
kxline=[0,1]; kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
kpoints=[Kx(:),Ky(:),Kz(:)];
critial=10^-5;
U=Us(i);
xinitial=[1,0,0,0;...
          0,0,0,0;...
          0,0,0,0;...
          0,0,0,1];
u1=1e-10;u2=0.5;
[xinitial,ni,si,efermi]=runhartreev5(g,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2);
g.onsite_modify(xinitial);
mz=[mz,abs((ni(1,1)-ni(1,2))-(ni(2,1)-ni(2,2)))]; %% Neel vectors

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"3DWeyl-bulk",0)
hold on;
plot(kpath,Energy(1,:)-efermi,'r-','LineWidth',2);
plot(kpath,Energy(2,:)-efermi,'b--','LineWidth',2);
plot(kpath,Energy(3,:)-efermi,'r-','LineWidth',2);
plot(kpath,Energy(4,:)-efermi,'b--','LineWidth',2);
end

%% how the magnetization evolves in an interacting Graphene from a AFM guess
Us=0:1:10;
lz=[];
egap=[];
for i=1:size(Us,2)
g.ham=g.iniham; 
stepmax=1000;
minstepmax=11;
step=1;
knum=101;
Electric_field_in_evpA=0;
kxline=[0,1]; kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
kpoints=[Kx(:),Ky(:),Kz(:)];
critial=10^-5;
U=Us(i);
xinitial=[1,0,0,0;...
          0,0,0,0;...
          0,0,0,0;...
          0,0,0,1];
u1=1e-10;u2=0.5;
[xinitial,ni,si,efermi]=runhartreev5(g,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2);
g.onsite_modify(xinitial);
lz=[lz,abs((ni(1,1)-ni(2,1))-(ni(1,2)-ni(2,2)))]; %% Neel vectors
[Energy,~,~]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
egap=[egap,min(Energy(3,:))-max(Energy(2,:))];
end
%
figure;
subplot(1,2,1);
plot(Us, egap, 'r-o');
% ylim([-0.2, 1.5]);
xlabel('U');
ylabel('E_{gap} (eV)');
title('Energy gap');
grid on; box on;
subplot(1,2,2);
plot(Us, lz, 'b-o');
ylim([-0.2, 2]);
xlabel('U');
ylabel('<lz>');
title('Magnetization vs Interaction');
grid on; box on;

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Run Hartree Fock for FM around VHs                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shift the Fermi Level to VHs
g.ham=g.iniham;
efermi=-1;
[Energy,kpath,kindex]= MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"3DWeyl-bulk",0)
ylim([-3,5])
%% Get the FM states by MF
U=4;
stepmax=1000;
minstepmax=11;
step=1;
knum=101;
Electric_field_in_evpA=0;
kxline=[0,1]; kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
kpoints=[Kx(:),Ky(:),Kz(:)];
critial=10^-5;
xinitial=[1,0,0,0;...
          0,0,0,0;...
          0,0,1,0;...
          0,0,0,0];
u1=1e-10;u2=0.375;
[xinitial,ni,si,efermi]=runhartreev5(g,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2);
g.onsite_modify(xinitial);
[Enk,Unk,kpath,kindex]=MTB.ham.get_bulk_bands_full(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Enk,nbands,efermi,kpath,labels,kindex,"3DWeyl-bulk",0)
hold on;
[sz_band,~]=get_band_sz(Unk,Enk);
for i=1:nbands
    scatter(kpath, Enk(i,:)-efermi, [], real(sz_band(i,:)), 'filled');
end
xlabel('Momentum'); ylabel('Energy'); set(gca,'XTick',[]);
xlim([min(kpath) max(kpath)]); box on; colorbar;colormap(blueWhiteRed(256));
clim([-max(abs(sz_band(:))) max(abs(sz_band(:)))]);
% colorbar 标签
c = colorbar;
c.Label.String = '<S_z>';
% c.Label.Interpreter = 'latex';
ylim([-3,5])

%% how the magnetization evolves in an interacting Graphene of a FM guess
Us=0:0.5:4;
mz=[];
for i=1:size(Us,2)
g.ham=g.iniham; 
stepmax=1000;
minstepmax=11;
step=1;
knum=101;
Electric_field_in_evpA=0;
kxline=[0,1]; kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
kpoints=[Kx(:),Ky(:),Kz(:)];
critial=10^-5;
U=Us(i);
xinitial=[1,0,0,0;...
          0,0,0,0;...
          0,0,1,0;...
          0,0,0,0];
u1=1e-10;u2=0.375;
[xinitial,ni,si,efermi]=runhartreev5(g,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2);
g.onsite_modify(xinitial);
mz=[mz,abs((ni(1,1)-ni(2,1))+(ni(1,1)-ni(2,2)))]; %% Neel vectors
% [Enk,Unk,kpath,kindex]=MTB.ham.get_bulk_bands_full(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% MTB.plot.plot_bands(Enk,nbands,efermi,kpath,labels,kindex,"3DWeyl-bulk",0)
end
figure;
plot(Us, mz, 'r-o');
% ylim([-0.2, 1.5]);
xlabel('U');
ylabel('<mz>');
title('Magnetization vs Interaction');   % ← 标题
grid on; box on;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Function Part                              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = blueWhiteRed(m)
%BLUEWHITERED Generate a blue–white–red diverging colormap.
%   cmap = BLUEWHITERED(m) returns an m-by-3 colormap that transitions
%   from blue (negative) to white (zero) to red (positive).
%   If m is omitted, defaults to 256.
    if nargin < 1, m = 256; end
    r = [(0:1/(m/2-1):1)  ones(1,m/2)];
    g = [(0:1/(m/2-1):1)  (1:-1/(m/2-1):0)];
    % Slightly dim the blue side by dividing by 1.1 for a softer low end
    b = [ones(1,m/2)/1.1  (1:-1/(m/2-1):0)];
    cmap = [b(:) g(:) r(:)];
end


function [sz_band, Enk] = get_band_sz(Unk, Enk)
%GET_BAND_SZ Compute <S_z> for each band at each k-point.
%   [sz_band, Enk] = GET_BAND_SZ(Unk, Enk)
%   Inputs
%     Unk : (dim x dim x Nk) eigenvector matrices at each k (columns are eigenstates)
%     Enk : (unused here except returned) energies; kept for API compatibility
%   Outputs
%     sz_band : (dim x Nk) expectation values <S_z> per band per k
%     Enk     : passed-through input (unchanged)
%
%   Note: Assumes spinful basis with ordering [orbital1(↑,↓), orbital2(↑,↓), ...].

    kpoints = size(Unk,3);
    dim = size(Unk,1);

    paulix = [0,1;1,0]; %#ok<NASGU>  % Pauli σ_x (not used; kept for reference)
    pauliy = [0,-1i;1i,0]; %#ok<NASGU>% Pauli σ_y (not used; kept for reference)
    pauliz = [1,0;0,-1];             % Pauli σ_z

    sz_band = zeros(dim, kpoints);
    sz = kron(eye(dim/2), pauliz);   % S_z operator in this basis (±1 eigenvalues)

    for k_idx = 1:kpoints
        psik = Unk(:,:,k_idx);               % eigenvectors at k (dim x dim)
        sz_val = diag(psik' * sz * psik);    % diagonal = <ψ_n|S_z|ψ_n>
        % Example of sorting two specific bands (commented out):
        % [sz_val(2:3),idx] = sort(sz_val(2:3), 'descend');
        sz_band(:,k_idx) = sz_val;
        % Example of reordering energies accordingly (commented out):
        % Enk(2:3,k_idx) = Enk(idx+1,k_idx);
    end
end


function result_matrix = find_neighbor_data(coords_cartesian, lattice_vectors, n, dimensionality)
%FIND_NEIGHBOR_DATA Find n-th neighbor pairs and return fractional offsets and distance.
%   result_matrix = FIND_NEIGHBOR_DATA(coords_cartesian, lattice_vectors, n, dimensionality)
%   Inputs
%     coords_cartesian : (N_atoms x 3) Cartesian coordinates within the unit cell
%     lattice_vectors  : (3 x 3) lattice basis vectors as rows or columns (consistent with multiplication)
%     n                : which neighbor shell to extract (1 = nearest, 2 = next nearest, ...)
%     dimensionality   : 2 for 2D lattice (ignore z separations), 3 for full 3D
%   Output
%     result_matrix    : (num_pairs x 6) matrix with rows [i, j, frax, fray, fraz, distance],
%                        where (frax,fray,fraz) are fractional offsets from atom i to the image of atom j.

    % 1) Convert Cartesian coordinates to fractional coordinates
    inv_lattice = inv(lattice_vectors);
    coords_fractional = coords_cartesian * inv_lattice;

    % 2) Build a supercell to capture neighbors across periodic images
    if dimensionality == 2
        % Extend only x and y directions
        [super_coords, super_indices] = construct_supercell_2d(coords_fractional, lattice_vectors, 1);
    elseif dimensionality == 3
        % Extend x, y, z directions
        [super_coords, super_indices] = construct_supercell(coords_fractional, lattice_vectors, 1);
    else
        error('Dimensionality must be 2 or 3');
    end

    % 3) Compute pairwise distances without minimum-image PBC correction (we already expanded)
    if dimensionality == 2
        % Consider only x-y planar distances
        dist_matrix_non_pbc = compute_distance_matrix_2d(super_coords, lattice_vectors, false);
    elseif dimensionality == 3
        % Full 3D distances
        dist_matrix_non_pbc = compute_distance_matrix(super_coords, lattice_vectors, false);
    end

    % 4) Extract unique distances (upper triangle), keep nonzero, group by rounding (tolerance bucket)
    raw_distances = triu(dist_matrix_non_pbc);
    % Optionally remove exact zeros (same atom); upper-tri already excludes i==j
    rounded_distances = round(raw_distances, 2);  % bucket distances to 2 decimals
    unique_distances = unique(rounded_distances, 'sorted');

    % 5) Select the n-th unique distance
    if n > length(unique_distances)
        error('n = %d exceeds the maximum neighbor shell.', n);
    end
    nth_distance = unique_distances(n);

    % 6) Find all pairs at that distance (within tolerance)
    [pair_i, pair_j] = find(abs(dist_matrix_non_pbc - nth_distance) < 1e-2);

    % 7) Build the result matrix
    result_matrix = find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors);
end


function unique_pairs = find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors)
%FIND_UNIQUE_NTH_NEIGHBORS Build unique n-th neighbor pairs with fractional offsets and distance.
%   unique_pairs = FIND_UNIQUE_NTH_NEIGHBORS(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors)
%   Input
%     pair_i, pair_j   : indices in the supercell arrays that form n-th neighbor pairs
%     super_coords     : (N_super x 3) fractional coords in the supercell
%     super_indices    : (N_super x 4) [atom_id_in_cell, i_shift, j_shift, k_shift]
%     nth_distance     : the selected neighbor distance (approximate)
%     lattice_vectors  : (3 x 3) lattice basis
%   Output
%     unique_pairs     : (num_pairs x 6) rows [i, j, frax, fray, fraz, distance]
%
%   Duplicates like (i->j) and (j->i) are removed.

    num_pairs = length(pair_i);
    all_pairs = zeros(num_pairs, 6);  % [i, j, frax, fray, fraz, dis]

    for k = 1:num_pairs
        i = pair_i(k); % index in supercell list
        j = pair_j(k); % index in supercell list

        % Map back to atom indices inside the home unit cell
        atom_i = super_indices(i, 1);
        atom_j = super_indices(j, 1);

        % Fractional displacement from i to j (includes supercell shift)
        delta_r = super_coords(i, :) - super_coords(j, :);

        % Convert to Cartesian to compute the actual distance
        delta_cartesian = delta_r * lattice_vectors;
        distance = sqrt(sum(delta_cartesian.^2));

        % Keep pairs close to the selected shell distance
        if abs(distance - nth_distance) < 1e-1
            % Store [i_cell, j_cell, fractional offset (j - i), distance]
            all_pairs(k, :) = [atom_i, atom_j, super_indices(j,2:end)-super_indices(i,2:end), nth_distance];
        end
    end

    % Remove duplicate rows (e.g., i->j and j->i)
    [~, unique_rows] = unique(all_pairs, 'rows');
    unique_pairs = all_pairs(unique_rows, :);
end


function dist_matrix = compute_distance_matrix(coords, lattice_vectors, pbc)
%COMPUTE_DISTANCE_MATRIX Pairwise distances in 3D (optionally with minimum-image PBC).
%   dist_matrix = COMPUTE_DISTANCE_MATRIX(coords, lattice_vectors, pbc)
%   Inputs
%     coords          : (N x 3) fractional coordinates
%     lattice_vectors : (3 x 3) lattice basis
%     pbc             : logical flag; if true, apply minimum-image convention
%   Output
%     dist_matrix     : (N x N) symmetric distance matrix (in Cartesian units)

    n = size(coords, 1);
    dist_matrix = zeros(n, n);

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                % Wrap to nearest image in fractional coordinates
                delta_r = delta_r - round(delta_r);
            end
            delta_cartesian = delta_r * lattice_vectors;
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % symmetry
        end
    end
end


function dist_matrix = compute_distance_matrix_2d(coords, lattice_vectors, pbc)
%COMPUTE_DISTANCE_MATRIX_2D Pairwise distances in 2D (ignore z component).
%   dist_matrix = COMPUTE_DISTANCE_MATRIX_2D(coords, lattice_vectors, pbc)
%   Inputs
%     coords          : (N x 3) fractional coordinates
%     lattice_vectors : (3 x 3) lattice basis
%     pbc             : logical flag; if true, apply minimum-image in 2D
%   Output
%     dist_matrix     : (N x N) distances in the x-y plane (Cartesian units)

    n = size(coords, 1);
    dist_matrix = zeros(n, n);

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r);
            end
            delta_cartesian = delta_r * lattice_vectors;
            delta_cartesian = delta_cartesian(:, 1:2); % keep x and y only
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j);
        end
    end
end


function [super_coords, super_indices] = construct_supercell(coords, lattice_vectors, scale)
%CONSTRUCT_SUPERCELL Build a 3D supercell by integer shifts in {-scale..scale}^3.
%   [super_coords, super_indices] = CONSTRUCT_SUPERCELL(coords, lattice_vectors, scale)
%   Inputs
%     coords          : (N_atoms x 3) fractional coordinates in the unit cell
%     lattice_vectors : (3 x 3) lattice basis (unused directly but kept for consistency)
%     scale           : number of copies in each ± direction
%   Outputs
%     super_coords    : ((2*scale+1)^3 * N_atoms) x 3 fractional coordinates
%     super_indices   : same rows, [atom_id_in_cell, i_shift, j_shift, k_shift]

    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            for k = -scale:scale
                offset = [i, j, k];
                super_coords = [super_coords; coords + offset]; %#ok<AGROW>
                for a = 1:n_atoms
                    super_indices = [super_indices; a, i, j, k]; %#ok<AGROW>
                end
            end
        end
    end
end


function [super_coords, super_indices] = construct_supercell_2d(coords, lattice_vectors, scale)
%CONSTRUCT_SUPERCELL_2D Build a 2D supercell (extend x and y only; z-shift = 0).
%   [super_coords, super_indices] = CONSTRUCT_SUPERCELL_2D(coords, lattice_vectors, scale)
%   Inputs/Outputs are analogous to CONSTRUCT_SUPERCELL, but k-shift is 0.

    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            offset = [i, j, 0]; % no extension along z
            super_coords = [super_coords; coords + offset]; %#ok<AGROW>
            for a = 1:n_atoms
                super_indices = [super_indices; a, i, j, 0]; %#ok<AGROW>
            end
        end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                    Hartree–Fock Driver Functions                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xinitial, ni, si, efermi] = runhartreev5(gs, knum, Kx, Ky, Kz, Electric_field_in_evpA, xinitial, stepmax, critial, U, u1, u2)
%RUNHARTREEV5 One outer HF cycle: relax onsite via quasi-Newton and compute observables.
%   [xinitial, ni, si, efermi] = RUNHARTREEV5(gs, knum, Kx, Ky, Kz, Efield, xinitial, stepmax, crit, U, u1, u2)
%   Inputs
%     gs         : structure with Hamiltonian and methods (expects gs.ham, gs.onsite_modify, etc.)
%     knum       : linear k-grid size (total Nk = knum^2)
%     Kx,Ky,Kz   : k-path components / grid definitions (passed to MTB.ham.get_bulk_plane_bands_add_electric)
%     Efield     : electric field in eV/Å
%     xinitial   : initial onsite matrix (vectorized dim x dim)
%     stepmax    : (unused here; kept for interface compatibility)
%     critial    : (unused here; kept for interface compatibility)
%     U          : Hubbard U parameter
%     u1,u2      : filling fractions defining an energy window for correlations
%   Outputs
%     xinitial   : relaxed onsite matrix (vectorized)
%     ni         : onsite densities (2 x Norb), row1 up, row2 down
%     si         : onsite spin vector components [Sx; Sy; Sz] (3 x Norb)
%     efermi     : Fermi level (by mid-gap between occupied and next state)
%
%   Note: This function relies on external MTB.ham methods for band calculation.

    nbands = size(gs.ham,1);

    xinitial = xinitial(:);
    objective = @(xinitial) one_step_hf_v3(gs, knum, Kx, Ky, Kz, Electric_field_in_evpA, xinitial, U, u1, u2);
    xinitial = quasi_newton(objective, xinitial);     % relax onsite via quasi-Newton

    % Apply relaxed onsite to the Hamiltonian
    xinitial = reshape(xinitial, nbands, nbands);
    gs.onsite_modify(xinitial);

    % Compute bands with electric field
    [Unk, Enk] = MTB.ham.get_bulk_plane_bands_add_electric(gs, Electric_field_in_evpA, Kx, Ky, Kz);
    Unk = reshape(Unk, [nbands, nbands, knum^2]);
    Enk = reshape(Enk, [knum^2, nbands]);

    % Total energy / occupied indices (by filling u2)
    [~, kindex, bandindex, efermi] = Total_energy(Enk, u2);

    % Onsite densities and spin expectation per site (from occupied states)
    [ni, si] = calonsite(Unk, kindex, bandindex);

    % Recompute Fermi energy as average of the two levels near filling
    efermi = calculate_ef_2(Enk, u2);
end


function xnew = one_step_hf_v3(gs, knum, Kx, Ky, Kz, Electric_field_in_evpA, xinitial_0, U, u1, u2)
%ONE_STEP_HF_V3 One objective evaluation for the quasi-Newton update.
%   Returns the updated onsite potential vector xnew given xinitial_0.
%   The update is a linear mixing of newly estimated onsite_U and current xinitial_0.

    nbands = size(gs.ham,1);

    xinitial_0 = reshape(xinitial_0, nbands, nbands);
    gs.onsite_modify(xinitial_0);

    [Unk, Enk] = MTB.ham.get_bulk_plane_bands_add_electric(gs, Electric_field_in_evpA, Kx, Ky, Kz);
    Unk = reshape(Unk, [nbands, nbands, knum^2]);
    Enk = reshape(Enk, [knum^2, nbands]);

    % Estimate correlation matrix in an energy window around Fermi set by (u1,u2)
    correlation = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);

    % Build onsite U matrix from correlations
    onsite_U = calculate_onsite_U(correlation, U);

    % Mix new onsite with previous (simple linear mixing)
    xnew = onsite_U * 0.8 + 0.2 * xinitial_0;
    xnew = xnew(:);
end


function C_total_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2)
%CALCULATE_CORRELATION_WITH_RANGE_AVG Correlation matrix averaged over k in [E(u1), E(u2)].
%   C_total_avg = CALCULATE_CORRELATION_WITH_RANGE_AVG(Enk, Unk, u1, u2)
%   - Finds energies E_low,E_high corresponding to fillings u1 and u2
%   - For each k, constructs a projector onto states with E in [E_low,E_high]
%   - Accumulates C(k) = conj(U_k) * W_k * U_k^T over k and averages

    [knum2, nbands] = size(Enk);
    C_total = zeros(nbands, nbands);

    % 1) Energy window by fillings
    E_low  = calculate_ef(Enk, u1);
    E_high = calculate_ef(Enk, u2);

    % 2) Accumulate over all k-points
    for m = 1:knum2
        E_k = Enk(m, :);         % (1 x nbands)
        U_k = Unk(:, :, m);      % (nbands x nbands)

        % Projector onto bands with energies in [E_low, E_high]
        W_k = diag(double(E_k >= E_low & E_k <= E_high));

        % Correlation contribution at this k
        C_m = conj(U_k) * W_k * transpose(U_k);

        C_total = C_total + C_m;
    end

    % 3) Average over k
    C_total_avg = C_total / knum2;
end


function onsite_U_old = calculate_onsite_U(correlation, U)
%CALCULATE_ONSITE_U Build onsite interaction matrix from correlation via a basis transform.
%   onsite_U_old = CALCULATE_ONSITE_U(correlation, U)
%   Inputs
%     correlation : (2*nbands x 2*nbands) correlation matrix in the spinful basis
%     U           : Hubbard U parameter
%   Output
%     onsite_U_old: onsite matrix (same size) after a block transform and masking
%
%   Notes:
%     - T = [0 -1; 1 0] is used as a 2x2 transform; T_large = kron(I_nbands, T)
%     - A mask of ones(2) on each orbital is applied afterwards
%     - The function returns the full matrix (diagonal extraction is commented out)

    total_bands = size(correlation, 1);
    if mod(total_bands, 2) ~= 0
        error('Correlation matrix must have even dimension (spin up/down per orbital).');
    end
    nbands = total_bands / 2;

    % Small transform matrix and its block-diagonal extension
    T = [0, -1; 1, 0];
    I_nbands = eye(nbands);
    T_large = kron(I_nbands, T);

    % Apply transform and mask
    mask_matrix = kron(I_nbands, ones(2));
    onsite_U_old = U * (T_large * correlation * T_large') .* mask_matrix;

    % Optionally reduce to diagonal only:
    % onsite_U_old = diag(diag(onsite_U_old));
end


function efermi = calculate_ef(Enk, u)
%CALCULATE_EF Estimate Fermi level corresponding to filling u (0..1).
%   efermi = CALCULATE_EF(Enk, u)
%   - Flattens energies and takes the maximum among the lowest ceil(u*N) values.

    total_states   = numel(Enk);
    occupied_states = ceil(total_states * u);
    efermi = max(mink(Enk(:), occupied_states));
end


function efermi = calculate_ef_2(Enk, u)
%CALCULATE_EF_2 Estimate Fermi level as midpoint between last occupied and first unoccupied.
%   efermi = CALCULATE_EF_2(Enk, u)

    total_states    = numel(Enk);
    occupied_states = ceil(total_states * u);
    efermi = (max(mink(Enk(:), occupied_states)) + ...
              max(mink(Enk(:), occupied_states + 1))) / 2;
end


function [ni, si] = calonsite(Unk, kindex, bandindex)
%CALONSITE Compute onsite densities and spin components from occupied states.
%   [ni, si] = CALONSITE(Unk, kindex, bandindex)
%   Inputs
%     Unk       : (Dim x Dim x Nk) eigenvectors; Dim = 2 * Nsites (spinful)
%     kindex    : indices of k points of occupied states
%     bandindex : indices of band numbers corresponding to kindex
%   Outputs
%     ni : (2 x Nsites) onsite density; row 1 = spin up, row 2 = spin down
%     si : (3 x Nsites) onsite spin components [Sx; Sy; Sz]
%
%   Notes:
%     - Uses parfor over occupied states to accumulate contributions.

    nki = zeros(size(Unk,1), 1);
    sitenum = size(Unk,2) / 2;

    sxki = zeros(sitenum, sitenum);
    syki = zeros(sitenum, sitenum);
    szki = zeros(sitenum, sitenum);

    paulix = [0,1;1,0];
    pauliy = [0,-1i;1i,0];
    pauliz = [1,0;0,-1];

    parfor i = 1:size(kindex,1)
        % Density contribution from occupied eigenstate |ψ>
        nki = abs(Unk(:, bandindex(i), kindex(i))).^2 + nki;

        % Reshape to (spin x site) for local spin expectation
        psik = reshape(Unk(:, bandindex(i), kindex(i)), [2, sitenum]);

        sxki = psik' * paulix * psik ./ 2 + sxki;  % <Sx> = σ_x/2
        syki = psik' * pauliy * psik ./ 2 + syki;  % <Sy> = σ_y/2
        szki = psik' * pauliz * psik ./ 2 + szki;  % <Sz> = σ_z/2
    end

    knum = size(Unk,3);          % total k-points
    ni = nki / knum;             % average over k
    ni = reshape(ni, [2, size(ni,1)/2]);  % first row up, second row down

    % Extract onsite spin components (diagonal in site space), average over k
    si = real([diag(sxki).'; diag(syki).'; diag(szki).'] ./ knum);
end


function [Etot, kindex, bandindex, efermi] = Total_energy(Enk, u)
%TOTAL_ENERGY Sum of occupied single-particle energies and occupied indices.
%   [Etot, kindex, bandindex, efermi] = TOTAL_ENERGY(Enk, u)
%   Inputs
%     Enk : (Nk x Nband) energies
%     u   : filling fraction (0..1)
%   Outputs
%     Etot      : total band energy per k (sum of occupied / Nk)
%     kindex    : k indices (row indices) of occupied states
%     bandindex : band indices (column indices) of occupied states
%     efermi    : Fermi energy estimated as max among occupied states

    [knum, ~] = size(Enk);

    % Collect lowest ceil(u * total_states) energies/indices
    [a, b] = mink(Enk(:), ceil(numel(Enk) * u));
    efermi = max(a, [], "all");

    % Map linear indices back to (row=k, col=band)
    row = mod(b, knum); row(row == 0) = knum;
    col = ceil(b ./ knum);

    Etot = sum(a, 'all') / knum;
    kindex = row;
    bandindex = col;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Quasi-Newton (Anderson-like) Accelerator             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xn = quasi_newton(func, x0, q, tol, history)
%QUASI_NEWTON Simple quasi-Newton/Anderson acceleration for fixed-point problems.
%   xn = QUASI_NEWTON(func, x0, q, tol, history)
%   Inputs
%     func    : function handle returning f(x) (same length as x)
%     x0      : initial guess (vector)
%     q       : history depth (default 20)
%     tol     : convergence threshold on max|Δx| (default 1e-5)
%     history : cell array to store iterates (optional)
%   Output
%     xn      : converged vector (or last iterate if not converged)
%
%   Notes:
%     - Performs a warm-up to build history matrices U and V.
%     - After each step, enforces nonnegative diagonal when reshaped as a square matrix
%       (model-specific constraint kept from original code).
%     - Adds small regularization if C is ill-conditioned.

    if nargin < 3
        q = 20;
    end
    if nargin < 4
        tol = 1e-5;
    end
    if nargin < 5
        history = {};
    end

    f0 = func(x0);
    n = length(f0);

    alpha = 1e-3;       % simple damping for warm-up
    U = zeros(n, q);
    V = zeros(n, q);

    xn = x0;
    fn = f0;

    disp('Starting quasi-Newton...');

    % ---- Warm-up: build U and V ----
    for i = 1:q
        fn = func(xn);
        xn_1 = xn - alpha * fn;

        % If x represents a square matrix, project diagonal to be nonnegative
        N = sqrt(length(xn_1));
        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        U(:, i) = fn - xn;
        V(:, i) = func(fn) - fn;
        xn = xn_1;
    end

    disp('Warmed up.');

    % ---- Main loop ----
    for i = 1:1000
        fn = func(xn);

        % Build C and RHS b
        C = transpose(U) * U - transpose(U) * V;
        b = transpose(U) * (xn - fn);

        % Regularize if singular/ill-conditioned
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        delta = V * (C \ b);
        xn_1 = fn - delta;

        % Project diagonal nonnegative when reshaped square (problem-specific)
        N = sqrt(length(xn_1));
        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        history{end+1} = xn_1; %#ok<AGROW>

        fn_1 = func(xn_1);
        U = [U(:, 2:end), fn_1 - xn_1];         %#ok<AGROW>
        V = [V(:, 2:end), func(fn_1) - fn_1];   %#ok<AGROW>

        dx = max(abs(xn_1 - xn));
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Finished at iteration %d\n', i);
            return;
        end

        xn = xn_1;
    end

    disp('DIDN''T CONVERGE!!!!');
end

