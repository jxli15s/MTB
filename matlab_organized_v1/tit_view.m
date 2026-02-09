% profile on
[xinitial, metaData] = flattenNestedCell(xinitial);        % 展平操作
xinitial = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, V, V, u1, u2, pairsU0, pairsU, pairsV);
xinitial = restoreNestedCell(xinitial, metaData);  % 复原操作
modifyHam(gs, xinitial, V, V, pairsU, pairsV)  %修改Ham
nbands=size(gs.ham,1);