function [vertices_ReT, faces_ReT, n_rem, ubelongcd, nfig, xdelta] = ...
    Re_tiling(vertices, faces, nCand, k_level, nfig, detail_polt)
%%
pre_compute_ReT;

%% 1.generateVertices and constrain points.
disp('Generate random Candidate Vertices on surface...')
[nameF_cand, vertices_cand, CSP_idx] = generateVertices2(mymesh, nCand, k_level);

% CSP
vertices_CSP = vertices(CSP_idx, :); n_CSP = length(CSP_idx);

% plot v_cand and CSP
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
hold on
plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'r*');axis equal
hold off; title('Candidate Points and CSP')
end

%% prapare for loop
vertices_cand = [vertices_cand; vertices_CSP];
nameF_cand = [nameF_cand, zeros(1, n_CSP)];

% compute the distance of faces to faces and faces to CSP
disp('computing FSD and VSD...')
[FSD, num_FSD] = FacesShortDistance(vertices, faces, 2*radius, hedge_face);
disp(['the minius number of num_FSD:', num2str(min(num_FSD))]);
[VSD, num_VSD] = VerticesShortDistance(vertices, faces, vertices_CSP, 2*radius);
disp(['the minius number of num_vSD:', num2str(min(num_VSD))]);

%% 2.loop k = 100
disp('Move candidate points on mesh, looping...')
for n_move = 1:20
    % (1).repulsive force
    forces = ComputeRepulsiveForce_adj(mymesh, nameF_cand, vertices_cand, FSD, VSD);
    nameF_cand3 = nameF_cand; % for forces test(test 1)

    % (2).move candidate points
    % we neednt consider the CSP when moving on mesh
    [vertices_cand2, nameF_cand2] = moveOnMesh_PLUS(mymesh,  vertices_cand(1:nCand, :), ...
        nameF_cand(1:nCand), forces(1:nCand, :));
    vertices_cand(1:nCand, :) = vertices_cand2;
    nameF_cand(1:nCand) = nameF_cand2;
end

% force_test; % test 1, to check the force perpendicular to the face normal
aplane = sum(forces(1:nCand, :) .* norm_face(nameF_cand3(1:nCand), :),2);
xdelta = max(abs(aplane));

% plot, after moving candidate points
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
hold on
plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'ro');axis equal
hold off
title('Candidate Points and CSP after moving')
end

%% UbelongCD
UbelongCD; % Re-Tiling to PPS

%% 3.doMutualTesselation
disp('Do MutualTesselation...')
faces_Mutual = doMutualTesselation2(mymesh, nameF_cand(1:nCand), vertices_cand);

vertices_Mutual = [vertices; vertices_cand(1:nCand,:)];
mesh_test; % test 2, check the candidate points in triangles or not.

% plot, MutualTesselation
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));
axis equal; title('Mutual Tesselation');
end
%% 4.RemovingOldVertices
disp('Removing Old Vertices...')
[vertices_ReT, faces_ReT, n_rem] = RemovingOldVertices_cpp2(...
    mymesh, vertices_cand(1:nCand, :), faces_Mutual, CSP_idx);

% plot 
if detail_polt
figure(nfig); nfig = nfig + 1;
trimesh(faces_ReT, vertices_ReT(:,1), vertices_ReT(:,2), vertices_ReT(:,3));axis equal;
title('mesh after kicking old points');
end

end