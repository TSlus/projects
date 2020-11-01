% ��ÿ�����Ӧ�����ҵ������ͳһ�� delaunayTriangulation
function [vertices_cand, nameF_cand] = moveOnMesh_PLUS(mymesh, vertices_cand, nameF_cand, forces) %#ok<INUSL>
mesh_to_VFdetail;

dam_ = 0.1;
% �����µ��λ��
vertices_cand_old = vertices_cand;
vertices_cand = vertices_cand + dam_ * forces; % P�� = P + kS

% �жϴ� push �ĵ��Ƿ���ԭ�������ڲ�
flag = isOnTriangle(nameF_cand, vertices_cand, vertices, faces, norm_face);  %#ok<USENS>
% 1.�����������ڲ��ĵ�ҲͶӰһ��
idx_not_out =  find(flag == 1);
for i = idx_not_out
    fi = nameF_cand(i);
    vertices_cand(i, :) = project_point_to_triangle(...
               vertices_cand(i, :), vertices(faces(fi,:),:), norm_face(fi,:));
end

% 2.����ȷ����push�ĵ�
idx_out =  find(flag == 0);% ��push ����������ĵ�Ķ�������
num_out = length(idx_out);
if num_out  == 0
    return;
end

%% Ѱ��ÿ��candidate �������λ��
for i = idx_out % i����candidate ������
    push_point = vertices_cand_old(i,:);
    push_point_towards = vertices_cand(i,:);
    
    for k = 1:20
        v1 = vertices(faces(nameF_cand(i), 1),:);
        v2 = vertices(faces(nameF_cand(i), 2),:);
        v3 = vertices(faces(nameF_cand(i), 3),:);
        vs = [v1; v2; v3]; vs2 = [v2; v3; v1];
        [t1, t2] = solve_two_cross_lines(vs, vs2, ...
            repmat(push_point, 3, 1), repmat(push_point_towards, 3, 1));
        t1_flag = (t1>0 & t1<1); t2_flag = (t2>0 & t2<1);
        t_idx = find(t1_flag & t2_flag); % �ҵ�һ�����������ཻ��
        if length(t_idx) ~= 1 % ����ֻ����0������ֹͣ push across
            vertices_cand(i, :) = project_point_to_triangle(...
                push_point, vertices(faces(nameF_cand(i),:),:), norm_face(nameF_cand(i),:));
            break;
        end
        t2 = t2(t_idx);
        cross_point = push_point + t2*(push_point_towards - push_point);
        
        % Ȼ������ת�ᣬ���� push_point
        e_num = t_idx;
        ek = [e_num, e_num+1];
        if e_num == 3
            ek = [3,1];
        end
        
        % ���� ek ȷ�� push ������һ����
        idx_vk = faces(nameF_cand(i),[ek(2), ek(1)]);
        fnum_changed = hedge_face(idx_vk(1), idx_vk(2));
        
        % �� ek_inv ��ת vout
        rotate_p = push_point_towards; % ��ת��
        rotate_line = vs(ek(2), :) - vs(ek(1), :);
        rotate_v0 = vs(ek(2), :);
        dot_ = dot(norm_face(nameF_cand(i),:), norm_face(fnum_changed,:));
        dot_(dot_ > 1) = 1; dot_(dot_ < -1) = -1;
        rotate_theta = acos(dot_);
        
        vp1 = ...
            point_rotate_line(rotate_p, rotate_line, rotate_v0, rotate_theta);
        vp2 = ...
            point_rotate_line(rotate_p, rotate_line, rotate_v0, - rotate_theta);
        vps = [vp1; vp2];
        
        dis_2 = dot(vps - rotate_v0, repmat(norm_face(fnum_changed,:),2,1), 2);
        dis_2 = abs(dis_2);
        idx_smaller = 1;
        if dis_2(1) > dis_2(2)
            idx_smaller = 2;
        end
        
        push_point_towards = vps(idx_smaller, :);
        push_point_old = push_point;
        push_point = cross_point + (1e-1)*(push_point_towards - cross_point);
        % push��ȥ�󣬼���Ƿ��� fnum_changed ���ڵ�ƽ��
        push_flag = isOnTriangle(fnum_changed, push_point, vertices, faces, norm_face);
        if ~push_flag
            vertices_cand(i, :) = project_point_to_triangle(...
                push_point_old, vertices(faces(nameF_cand(i),:),:), norm_face(nameF_cand(i),:));
            break;
        end
%         vertices_cand(i, :) = push_point;
        vertices_cand(i, :) = project_point_to_triangle(...
               push_point, vertices(faces(fnum_changed,:),:), norm_face(fnum_changed,:));
        nameF_cand(i) = fnum_changed;
    end
% ����candidate �������λ�õõ�����
end

end
