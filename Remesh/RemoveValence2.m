% 将remesh过后的网格去除度数是2的情况
function [vertices, faces] = RemoveValence2(vertices, faces)

nf = size(faces,1); 
faces_turn = faces(:, [2,3,1]);
f_name = [1:nf, 1:nf, 1:nf]';
half_face = sparse(faces(:), faces_turn(:), f_name);

np = size(vertices,1);
vertex_valence(np) = 0;
for P = 1:np
    neighbor_P = find(half_face(:,P));
    vertex_valence(P) = length(neighbor_P);   
end

P2 = find(vertex_valence == 2);

delet_f = zeros(nf, 1);
adj_idx = 1:np;
for i = P2
    [ri, ~] = find(faces == i);
    delet_f(ri) = 1;
    adj_idx((i+1):np) = adj_idx((i+1):np) - 1;
end
faces_new = faces(delet_f == 0,:);

faces = adj_idx(faces_new);
vertices = removerows(vertices, P2);
vertex_valence = removerows(vertex_valence', P2);

min_valen = min(vertex_valence);
disp(['The minimum degree:',num2str(min_valen)]);
end