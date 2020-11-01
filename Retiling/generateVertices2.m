function [nameF_cand, vertices_cand, CSP_idx] = generateVertices2(mymesh, nCand, k_level) 
mesh_to_VFdetail;

% add constrain points
% we set the point which has high curvature as CSP
disp('Add constrain points...')
[CSP_idx, abs_Cmean, abs_Cgaussian] = AddConstrainPoint(mymesh, k_level);
% 根据顶点的曲率计算面的曲率

abs_Cmean_f = sum(abs_Cmean(faces), 2);
abs_Cgaussian_f = sum(abs_Cgaussian(faces), 2);

% 选出面数0.2 的点作为高曲率面
n_face_HC = ceil( nf * 0.1);
[~,f_1] = maxk(abs_Cmean_f, n_face_HC);
[~,f_2] = maxk(abs_Cgaussian_f, n_face_HC);

f_HC = unique([f_1; f_2]); n_f_HC = length(f_HC);
% nCand 中有一部分点需要在高曲率出产生

% if nCand > 3 * n_f_HC
%     n_cand_HC = 3 * n_f_HC;
%     everyN = 3;
% elseif nCand > 2 * n_f_HC
if nCand > 2 * n_f_HC
    n_cand_HC = 2 * n_f_HC;
    everyN = 2;
elseif nCand > n_f_HC
    n_cand_HC = n_f_HC;
    everyN = 1;
else
    n_cand_HC = 0;
    everyN = 0;
end
%

face_highCurv = faces(f_HC,:);
vcenter_f_HC = (vertices(face_highCurv(:,1),:) + vertices(face_highCurv(:,2),:) + vertices(face_highCurv(:,3),:))/3;
all_HC_cand = (repmat(vcenter_f_HC, 3, 1) + vertices(face_highCurv(:),:))/2;
v_cand_HC = all_HC_cand(1:(n_f_HC * everyN), :);
temp = repmat(f_HC(1:n_f_HC), everyN, 1);
nameF_cand_HC = temp(:);

nCand = nCand - n_cand_HC;
% 这里的nCand已经除去两类点
vertices_cand = zeros(nCand, 3);
[rand_num, nameF_cand, s_sum_i] = GRPONN(si, nCand);

% 利用rand_num计算权重；
u1 = (rand_num - s_sum_i(nameF_cand)) ./ si(nameF_cand);
u2 = (1 - u1) .* rand(nCand, 1);
u3 = 1 - u1 - u2;

u = [u1, u2, u3];
if ~isempty(find(u<=0 | u>=1,1))
    disp('candidate points 权重u取值超出(0,1)，重新生成。')
end

% 生成随机点的坐标
for i = 1:nCand
    vertices_cand(i, :) = [u1(i), u2(i), u3(i)] * vertices(faces(nameF_cand(i),:),:);
end

nameF_cand = [nameF_cand_HC; nameF_cand]';
vertices_cand = [v_cand_HC; vertices_cand];

end