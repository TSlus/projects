%% ԭ�����LoopSurface����ͬ��ϸ�����κ�ģ�

% 1����������������LoopSurface��ԭ������Ҫϸ�����Σ�
% 2���ҵ�ԭ���������κ�ϸ�����κ�������֮��Ĺ�ϵ��hedge_face��
% 3���õ�ԭ���������κ�LoopSurface���ϵ��ԭ�����LoopSurface��
function loop_point = mesh_connect_LoopSurf(vertices, faces)
%% 1.LoopSurfaces
% load('model.mat');
faces_old = faces;
% ��������ϸ��
loopTimes = 2;
[vertices, faces, hedge_face] = myLoop(vertices, faces, loopTimes);
% trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal

% ����loopSurface
% ע�������loop_point��ϸ�ֺ������ζ�Ӧ
loop_point = LoopSurfaceSample(vertices, faces); 

%% 2.ͨ��ԭ���������κ�ϸ�����κ�������֮��Ĺ�ϵ
% �õ�ԭ�����LoopSurface
nf_old = size(faces_old, 1);
nf_new = size(faces, 1);

x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
R = [1:size(faces, 1), 1:size(faces, 1), 1:size(faces, 1)]';
he_f = sparse(X, Y, R);
connect{nf_old} = []; % nf ��ϸ��ǰfaces����
nEvery = 100; % loopϸ������loop_pointÿ�������βɵ����100��
% ���ѭ����������;����Ӧ����
fnew_to_fold = zeros(nf_old, 16);
tf = ones(nf_old,1);
for i = 1:nf_new
    fi = hedge_face(faces(i,1),faces(i,2));
    fnew_to_fold(fi, tf(fi)) = i;
    tf(fi) = tf(fi) + 1;
end
% �������ϵĵ㣬��Ӧ�������� 

for i = 1:nf_old
    cell_p = loop_point(fnew_to_fold(i,:));
    v_fi = zeros(16*nEvery, 3); tv = 1;
    for j = 1:16
        k = size(cell_p{j}, 1);
        v_fi(tv:tv+k-1, :) = cell_p{j};
        tv = tv + k;
    end
    v_fi =  v_fi(1:tv-1, :);
    connect{i} = v_fi;
end
loop_point = connect; 
%%
% for i=1:nf
%     [row, col] = find(hedge_face == i);
%     fi = he_f(sub2ind(size(he_f), row, col));
%     fi = unique(fi); n_fi = length(fi);
%     cell_p = loop_point(fi);
%     % ����ЩԪ���еĵ�ȡ����
%     v_fi = zeros(n_fi*nEvery, 3);
%     t = 1;
%     for j = 1:n_fi
%         k = size(cell_p{j}, 1);
%         v_fi(t:t+k-1, :) = cell_p{j};
%         t = t+k;
%     end
%      v_fi =  v_fi(1:t-1, :);
%      connect{i} = v_fi;
% end
% % ע�������loop_point��ϸ��ǰ�����ζ�Ӧ
% loop_point = connect; 

end
