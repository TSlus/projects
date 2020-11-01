function [vertices, faces, hedge_face] = myLoop(vertices, faces, loopTimes)
% ��������loopϸ��
%% ϸ�ֹ���
% �����loopϸ�ּ���һЩ����
% 1����ϸ�ֶ�κ�������棬��ԭ���������Ӧ���˴�û�У�����2���õ���
% 2��ÿ��ԭʼ��������Щϸ�ֵ�

% ����
% vertices : np x 3
% faces    : nf x 3

% �����
% vertices : np x 3��ϸ�ֽ���ĵ�
% faces    : nf x 3��ϸ�ֽ������
% hedge_face��loopϸ�ֺ�ģ����--�棩���ݽṹ��
%             �����ҵ�ÿ��ԭʼ��������Щloop�㣬��37-38��

%% mesh
nf = size(faces, 1);

%% face name
% ��ÿ��loop����������������Ϳ��ԴﵽĿ��1��
% ����Ŀ��1���ҵ�ͬһ��ԭʼ���loop�㣬�Ϳ���ʵ��Ŀ��2��
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; F = [1:nf, 1:nf, 1:nf]';
hedge_face_old = sparse(X, Y, F);
%%
% trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
for i = 1:loopTimes
    [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old);
    hedge_face_old = hedge_face;
end
% h = trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal
% set(h,'EdgeColor',[1.0, 0.66, 0.49],'LineWid',0.0001); 
% %% ���õõ���hedge_face���鿴ԭʼ���Ӧloop����Щ�µ�
% hold on;
% [row, col] = find(hedge_face == 3);
% points = unique([row; col]);
% hold on
% plot3(vertices(points,1),vertices(points,2),vertices(points,3),'r*');

end