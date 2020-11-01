% Ѱ�ҵ������㣬û�ж���
function [oneRingP, vertex_valence] = NearP(faces)

% faces -  nf * 3
% oneRingP - ���򶥵�
% vertex_valence - ����Ķ�

%���-��id
nf = size(faces,1);
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
f_name = [1:nf, 1:nf, 1:nf]';
half_face = sparse(X, Y, f_name);

np = max(max(faces));
oneRingP{np} = [];
vertex_valence(np) = 0;
for P = 1:np
    neighbor_P = find(half_face(:,P));%���ط���Ԫ��λ�ã���P������㣨û��˳��
    vertex_valence(P) = length(neighbor_P); 
    oneRingP{P} = neighbor_P; 
end
end
