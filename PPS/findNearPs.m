function [oneRingP, vertex_valence] = findNearPs(faces)

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
    %��ʼ������һ����
    neighbor_P = find(half_face(:,P));%���ط���Ԫ��λ�ã���P������㣨û��˳��
    vertex_valence(P) = length(neighbor_P);

    a = neighbor_P(1);
    oneRingP{P}(1) = a;
    t = 2;
    neighbor_res = neighbor_P(2:end);
    %�� b ����һ����
    while  ~isempty(neighbor_res)        
        a_nextP = find(half_face(a,neighbor_res));
        
        k = 1;
        b = neighbor_res(a_nextP(k)); %˳ʱ�뻹����ʱ�벻֪��
        while half_face(a,b) ~= half_face(P, a) % ������ʱ�����������Ҫ����
            k = k+1;
            b = neighbor_res(a_nextP(k));    % �������û��oneRing������ֱ����!!!
        end
        oneRingP{P}(t) = b; t = t+1;
        a = b; 
        neighbor_res = neighbor_res(neighbor_res ~= b);
    end
end
end
