addpath(genpath('Remesh')); 
addpath(genpath('Retiling')); addpath(genpath('PPS'))

%% the original mesh
nfig = 1;
% plot
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
title('the original mesh')

%% do remesh
disp('doing Remesh...')
[vertices, faces] = DoRemesh(vertices, faces);
% plot
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
title('the mesh after remesh')
end

%% do Re_Tiling
disp('======== Do Re-Tiling ========')
[vertices_ReT, faces_ReT, n_rem, ubelong, nfig, xdelta] = ...
    Re_tiling(vertices, faces, nCand, k_level, nfig, detail_polt);
% note：vertices_ReT中前n_rem个点是原始顶点

%% do PPS
% 1.loopSurface
disp('======== Do PPS ========')
disp('Constructing LoopSurface...')
loop_point = mesh_connect_LoopSurf(vertices, faces);

disp('Computing cadidate-PPS...')
vcdPPS = zeros(size(ubelong, 1), 3);

idx_odd = ubelong(:, 15) == 0; % 有部分点可能不被Omega区域覆盖
a = ubelong(idx_odd, 1); vap = vertices(a,:);
b = ubelong(idx_odd, 2); vbp = vertices(b,:);
c = ubelong(idx_odd, 3); vcp = vertices(c,:);
wights = ubelong(idx_odd, 4:6);
vcdPPS(idx_odd, :) = wights(:,1).*vap + wights(:,2).*vbp + wights(:,3).*vcp;

idx_ord = ubelong(:, 15) == 1;
ube_part = ubelong(idx_ord,:);

% 2.Retiling and PPS
vcdPPS(idx_ord, :) = surfaceConstruction(vertices, faces, loop_point, ube_part);
vertices_ReT(n_rem+1:end, :) = vcdPPS;
vertices_final = vertices_ReT; faces_final = faces_ReT;

% plot
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces_final, vertices_final(:,1), vertices_final(:,2), vertices_final(:,3));axis equal;
title('Retiling & PPS no remesh');
end

%% Remesh
disp('======== Remesh ========')
disp('Remeshing...')
[verticesr, facesr] = DoRemesh(vertices_final, faces_final);

% plot

figure(nfig); nfig = nfig + 1;
trimesh(facesr, verticesr(:,1), verticesr(:,2), verticesr(:,3));axis equal;
title('Retiling & PPS');

% cpp_result