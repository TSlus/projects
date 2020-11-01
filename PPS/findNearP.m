function oneRingP = findNearP(faces, P)

% faces -  numF * 3
% P - index

% ��� - ��
nf = size(faces,1);
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; 
f_name = [1:nf, 1:nf, 1:nf]';
half_face = sparse(X, Y, f_name);

%��ʼ������һ����
neighbor_P = find(half_face(:,P)); %���ط���Ԫ��λ�ã���P������㣨û��˳��
a = neighbor_P(1);
% �洢�����
oneRingP(length(neighbor_P)) = 0;
% length(neighbor_P)
oneRingP(1) = a;
t = 2;
neighbor_res = neighbor_P(2:end);

%����һ����
while  ~isempty(neighbor_res)
    a_nextP_idx = find(half_face(a, neighbor_res));
    k = 1; 
    b = neighbor_res(a_nextP_idx(k));
    while half_face(a,b) ~= half_face(P, a) 
        k = k + 1;
        b = neighbor_res(a_nextP_idx(k));    
    end
    oneRingP(t) = b;
    t = t + 1;
    a = b; neighbor_res = neighbor_res(neighbor_res ~= b);
end
end
