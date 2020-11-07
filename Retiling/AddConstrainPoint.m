function [CSP_idx,abs_Cmean, abs_Cgaussian] = AddConstrainPoint(mymesh, k_level)   %#ok<INUSL>
mesh_to_VFdetail; % ��mesh�е�������ɢ��

% �����ʽϴ�ĵ���ΪCSP
FV.vertices = vertices;
FV.faces = faces;
[Cmean, Cgaussian, ~,~,~,~] = GetCurvature(FV, true);
abs_Cmean = abs(Cmean);
abs_Cgaussian = abs(Cgaussian);

[~, CSP_idx1] = maxk(abs_Cmean, ceil(k_level/2));
[~, CSP_idx2] = maxk(abs_Cgaussian, ceil(k_level/2));

CSP_idx = sort(unique([CSP_idx1; CSP_idx2]));
end