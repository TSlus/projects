function [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old)
nf = size(faces, 1); np = size(vertices, 1);
nhe = 3 * nf; % �����
ne = nhe / 2;
if mod(nhe, 2) ~= 0
    warning('����������Ŀ����ż���������Ƿ������')
    return;
end

% ��߽ṹ
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
% ���
% he = zeros(nhe, 2);
% he(:,1) = X; he(:,2) = Y;
% ���--��
F = [1:nf, 1:nf, 1:nf]';
hedge_face = sparse(X, Y, F);
% ���--�Ե�
Z = [x3; x1; x2];
hedge_third_p = sparse(X, Y, Z);
% ���--�����ߵ�����
hedge_oppo_name = sparse(Y, X, (1:3*nf)');

% ÿ��ԭʼ����Ķ�
v_valence(np) = 0;
nearPs{np} = [];
for P = 1:np
    neighbor_P = find(hedge_face(:,P));
    nearPs{P} = neighbor_P;           % P�������
    v_valence(P) = length(neighbor_P);  % P�Ķ�
end

% ˼·�����ѭ������ÿ������ѳ��ĸ�
flag_he = zeros(nhe, 1); % ��¼������Ƿ�������µ㣬ͬʱ������¼����ĵ�
% hedge_midpoint = hedge_face; % ��ÿ������ϲ���һ����

% ��ѭ��
idxp = np + 1; % ������µ�����
faces_new = zeros(4*nf, 3);
vertices_new = zeros(np + ne, 3);
flag_v = vertices_new(:, 1); % ��¼�����Ƿ��Ѿ�����
face_name = zeros(12*nf, 3); % ��¼�²�����ߣ�����ԭʼ��
for i = 1:nf
    a = faces(i, 1); b = faces(i, 2); c = faces(i, 3);
    %% �����µ�k1, k2, k3
    % add_point;
    % ��һ��������
    if flag_he(i) % �����һ��������Ѿ������˵�
        k1 = flag_he(i);
    else          % �����һ������ϻ�û�����µ�
        k1 = idxp;
        flag_he(i) = idxp;
        % ������ҲҪ������Ӧ�ĵ�
        name = hedge_oppo_name(a, b);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    % �ڶ���������
    if flag_he(i + nf) % ����ڶ���������Ѿ������˵�
        k2 = flag_he(i + nf);
    else          % ����ڶ�������ϻ�û�����µ�
        k2 = idxp;
        flag_he(i + nf) = idxp;
        % ������ҲҪ������Ӧ�ĵ�
        name = hedge_oppo_name(b, c);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    % ������������
    if flag_he(i + 2*nf) % �����һ��������Ѿ������˵�
        k3 = flag_he(i + 2*nf);
    else          % �����һ������ϻ�û�����µ�
        k3 = idxp;
        flag_he(i + 2*nf) = idxp;
        % ������ҲҪ������Ӧ�ĵ�
        name = hedge_oppo_name(c, a);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    %%
    % 1.������
    face_add  = [a, k1, k3;
        k1, b, k2;
        k3, k2, c;
        k1, k2, k3];
    faces_new(4*i - 3:4*i,:) = face_add;
    
    y1 = [face_add(:,1),face_add(:,2)]; y2 = [face_add(:,2),face_add(:,3)];
    y3 = [face_add(:,3),face_add(:,1)];
    face_name((12*i-11):12*i, [1,2]) = [y1; y2; y3];
    face_name((12*i-11):12*i, 3) = hedge_face_old(a, b);
    
    %%
    % 2.���µ�
    % update_vertex;
    % a��
    if flag_v(a) == 0
        flag_v(a) = 1;
        va = vertices(a,:);
        u = v_valence(a);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(a,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{a},:);
    end
    
    % b��
    if flag_v(b) == 0
        flag_v(b) = 1;
        va = vertices(b,:);
        u = v_valence(b);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(b,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{b},:);
    end
    
    % c��
    if flag_v(c) == 0
        flag_v(c) = 1;
        va = vertices(c,:);
        u = v_valence(c);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(c,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{c},:);
    end
    
    %
    % k1
    if flag_v(k1) == 0
        flag_v(k1) = 1;
        v1 = vertices(a,:); v2 = vertices(b, :); v3 = vertices(c, :);
        d = hedge_third_p(b, a); v4 = vertices(d,:);
        vertices_new(k1, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    % k2
    if flag_v(k2) == 0
        flag_v(k2) = 1;
        v1 = vertices(b,:); v2 = vertices(c, :); v3 = vertices(a, :);
        d = hedge_third_p(c, b); v4 = vertices(d,:);
        vertices_new(k2, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    % k3
    if flag_v(k3) == 0
        flag_v(k3) = 1;
        v1 = vertices(c,:); v2 = vertices(a, :); v3 = vertices(b, :);
        d = hedge_third_p(a, c); v4 = vertices(d,:);
        vertices_new(k3, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    
end

vertices = vertices_new;
faces = faces_new;
hedge_face = sparse(face_name(:, 1), face_name(:, 2), face_name(:, 3));
end


