%��������������
% clear;clc;
% %% ���ݼ���
% load('model.mat')
% 
% %% ����Loop2_point
% load('model_loop2_A.mat');
function surfaceP = surfaceConstruction(vertices, faces, loop_point, ubelong)
% function surfaceP = surfaceConstruction(vertices, faces, loop_point)
%%
pre_compute_PPS; % Ԥ����

%% ֱ�Ӽ���������Ϣ
% ContrlPoint_Omega;
ctrPs = ContrlPoint_Omega2(vertices,vf_sparse,...
    v_valence, oneRingPs, loop_point, n_choose_k);
% ubelongGenerate;

%% 
%����ubelong�����жϳ�ÿ�������ĸ�Omegau��
nSample = size(ubelong,1);
surfaceP = zeros(nSample,3);
% load('n_choose_k.mat');
ti = 1;
% ���������㣨����ʱ�䣩
for kk = 1:nSample   
    Jumat = [ubelong(kk,[2, 11]);ubelong(kk,[3,14])];
    Ju = [ubelong(kk,1); Jumat(Jumat(:,2)==1, 1)];
    
    phi_vup_mat = [ubelong(kk,9:11); ubelong(kk,12:14)];
    phi_vup =  [ubelong(kk,7:8); phi_vup_mat(phi_vup_mat(:,3)==1,[1,2])];
    phi_vup = phi_vup';
    
    %����Ȩϵ��w(v,u)
    gammaup = zeros(1,length(Ju)); 
    tnorm = sum(abs(phi_vup).^2,1).^(1/2);
    H2 = cos(pi./v_valence(Ju)); H1 = 0.25 * H2;%ÿ��v��Ӧ��ͬ��H,tnorm
    gammaup(tnorm <= H1) = 1;
    gammaup(tnorm >= H2) = 0;
    H = (tnorm - H1)./(H2 - H1);
    perfidx = find(tnorm > H1 & tnorm < H2);
    
    ss = 1./sqrt(1-H(perfidx))-1./sqrt(H(perfidx));
    gammaup(perfidx)=1./(1 + exp(2*ss));
    wvu = gammaup / sum(gammaup);%�õ���λȨ��
    
    %��phi_vup�еĵ��������
    Lj = abs(cos(pi./v_valence(Ju)));
    phi_vup = (phi_vup + Lj) ./ (2 * Lj);%��һ��
    
    %����Ӧ��Ju�ϼ���Bezier��ֵ
    pi_v = zeros(length(Ju),3);
    for jj=1:length(Ju)
        mbctr = ctrPs{Ju(jj)};
        %Q��ÿһ������ά��
        m = v_valence(Ju(jj))+1;
        cB = zeros(1, m + 1);
        cB((0:m)+1) = n_choose_k(m+1, (0:m)+1);
        u = phi_vup(:,jj);
        Bu=cB.*(u.^(0:m)).*((1-u).^(m-(0:m)));
        Buv=Bu(1,:)'.* Bu(2,:);
        pi_v(jj,:) = Buv(:)' * mbctr;      
    end
    
    %�����ͼ����(u,p)
    theta_up = sum(wvu'.* pi_v,1);
    surfaceP(ti,:) = theta_up;
    ti = ti + 1;
end

% %% ��ͼ
% figure(1)
% plot3(surfaceP(:,1),surfaceP(:,2),surfaceP(:,3),".");axis equal

end
