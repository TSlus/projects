% ��������loopϸ��
%% ϸ�ֹ���
% �����loopϸ�ּ���һЩ����
% 1����ϸ�ֶ�κ�������棬��ԭ���������Ӧ
% 2��ÿ��ԭʼ��������Щϸ�ֵ�
%% ����
clear;clc;
load('model.mat'); faces_old = faces; vertices_old = vertices;
load('loop_point_model.mat') %���ڽ�loop�����Ӧ��ԭʼ����������

% load('bronze.mat'); faces_old = faces; vertices_old = vertices;
% load('loop_point_bronze.mat')

% vertices : np x 3
% faces    : nf x 3
%% mesh
nf = size(faces, 1); np = size(vertices, 1);

%% face name
% ��ÿ��loop����������������Ϳ��ԴﵽĿ��1��
% ����Ŀ��1���ҵ�ͬһ��ԭʼ���loop�㣬�Ϳ���ʵ��Ŀ��2��
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; F = [1:nf, 1:nf, 1:nf]';
hedge_face_old = sparse(X, Y, F);
%%
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
for i = 1:2
    [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old);
    hedge_face_old = hedge_face;
end
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal

%% 1.���õõ���hedge_face���鿴ԭʼ���Ӧloop2����Щ�µ�
hold on;
[row, col] = find(hedge_face == 3);
points = unique([row; col]);
hold on
plot3(vertices(points,1),vertices(points,2),vertices(points,3),'*');

%% 2.��ԭʼ��������ϸ�������ϲ�������ϵ����
% һ��ԭʼ�棬��Ӧϸ�����κ��16����

% ϸ�ֺ�İ��--��
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
R = [1:size(faces, 1), 1:size(faces, 1), 1:size(faces, 1)]';
he_f = sparse(X, Y, R);

connect{nf} = []; % nf ��ϸ��ǰfaces����
nEvery = 100; % loopϸ������loop_pointÿ�������βɵ����100��
for i=1:nf
    [row, col] = find(hedge_face == i);
    fi = he_f(sub2ind(size(he_f), row, col));
    fi = unique(fi); n_fi = length(fi);
    cell_p = loop_point(fi);
    % ����ЩԪ���еĵ�ȡ����
    v_fi = zeros(n_fi*nEvery, 3);
    t = 1;
    for j = 1:n_fi
        k = size(cell_p{j}, 1);
        v_fi(t:t+k-1, :) = cell_p{j};
        t = t+k;
    end
     v_fi =  v_fi(1:t-1, :);
     connect{i} = v_fi;
end
loop_point = connect;

%% ��ͼ��֤
figure(2)
faces = faces_old; vertices = vertices_old;
trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3));axis equal
hold on;
v1 = vertices(faces(1,:),:);
plot3(v1(:,1), v1(:,2), v1(:,3),'r*');
v2 = loop_point{1};
plot3(v2(:,1), v2(:,2), v2(:,3),'b.');
