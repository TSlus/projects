% Ԥ����
numP = size(vertices,1);numF = size(faces,1);
%% ������Ϣ
[oneRingPs, v_valence] = findNearPs(faces); %����Ԫ����ʽ�������
%% ����������Ƭ��������ڵ�һ����
tri_point_valence = v_valence(faces);
[irregular_f,extra_v] = find(tri_point_valence ~= 6);
irregular_f2 = irregular_f(extra_v == 2);
irregular_f3 = irregular_f(extra_v == 3);
faces(irregular_f2, :) = faces(irregular_f2, [2,3,1]); 
faces(irregular_f3, :) = faces(irregular_f3, [3,1,2]); 

% %% ����faces�����ж�
% % �ж��Ƿ�����ȵ�������ĵ�һ������
% flag1 = size(find(v_valence(faces(:,1)) ~= 6),2); %��һ�ж�������6����
% flag2 = find(v_valence(faces(:,2)) ~= 6 | v_valence(faces(:,3)) ~= 6);
% if (flag1~= length(irregular_f)) || size(flag2,2)
%     disp('�����ε�������')
% end

%% ���������������
regular_f = setdiff(1:numF, irregular_f);      %������������������
faces_regular = faces(regular_f,:);

% �����������棺1.������3���档2.��������3���档
tri_point_valence = v_valence(faces);
r3 = (tri_point_valence(irregular_f, 1) == 3);
irregular_f_3 = irregular_f(r3);
irregular_f = irregular_f(~r3);
faces_irregular = faces(irregular_f,:);

%% ��������������
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; R = [1:numF, 1:numF, 1:numF]';
vf_sparse = sparse(X, Y, R);

%% n_choose_k
max_val = max(v_valence);
n_choose_k = sparse(max_val + 2, max_val + 2);
for i = 0:max_val + 1
    for j = 0:max_val + 1
    if i >= j
        n_choose_k(i+1, j+1) = nchoosek(i,j);
    end
    end
end

%% �ȼ��� N=6����Ŀ��ƶ���
cr_idx = controlP_idx_regular_face(faces_regular, oneRingPs); % nf1 x 12
% �ټ���N!=6�Ŀ��ƶ���
cir_idx = controlP_idx_irregular_face(faces_irregular, oneRingPs, v_valence);

%% ϸ�־�������
myEigen = EIGENSTRUCT; % NMIN = 4; NMAX = 100;
