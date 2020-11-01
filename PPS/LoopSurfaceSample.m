% Loop�������
%%
% clear; clc;
% load('model_loop2_A.mat') %����ϸ�����εĽ��

function loop_point = LoopSurfaceSample(vertices, faces)
% ���룺һ������ϸ�����ε�����

pre_compute_loopSurf;
nEveryF = 10; % ÿ�����ϲ�15����
%% 1��������Ƭ���飺���������������Ρ�������������������
irregular_f_order = sort(irregular_f);       %����������������
is_isolate = irregular_f_order(1:end - 1) - irregular_f_order(2 : end);
if find(is_isolate == 0)
    disp('����㲻����');
end
n_irregular_f = length(irregular_f);
n_regular_f = length(regular_f);

%% 2����������û��������������
%ȷ��12���������������������ֱ�Ӽ���

%ȷ��12�������λ��
cr_idx = cr_idx';
C0_regular = vertices(cr_idx(:), :)'; % 3 x 12*rf

% ���ı�ʾ��ʽ
b_vw_m; % ��������
temp_mat = reshape(1:12*n_regular_f, 12, n_regular_f)';
C0_regular = C0_regular(:, temp_mat(:));
regular_patch = reshape(C0_regular, 3*n_regular_f, 12) * b_vw'; % C0 * b

%% 3���������������������
% ����һ���ṹ�壬��ϸ�־��� A �������ṹ
myEigen = EIGENSTRUCT; % NMIN = 4; NMAX = 100;

% ���ڴ��������������Σ��ҵ���Χ�Ŀ��Ƶ�C_0��(N+6) x 3��
% ���Ƶ�C_0�����ڲ�ͬ N���ҵ�(N+6)�����ƶ��㣬����ֵ��һ��Ԫ�� 1 x NMAX
cir_idx = controlP_idx_irregular_face(faces_irregular, oneRingPs, v_valence);

% ������������patch
irregular_surf = EvalSurf(vertices, faces, v_valence, irregular_f, cir_idx, myEigen, nEveryF);
loop_point_arran; % ����������� n x 3

end