function [CSP_idx,abs_Cmean, abs_Cgaussian] = AddConstrainPoint(mymesh, k_level)   %#ok<INUSL>
mesh_to_VFdetail; % 将mesh中的数据离散化

% 将曲率较大的点设为CSP
FV.vertices = vertices;
FV.faces = faces;
[Cmean, Cgaussian, ~,~,~,~] = GetCurvature(FV, true);
abs_Cmean = abs(Cmean);
abs_Cgaussian = abs(Cgaussian);

[~, CSP_idx] = maxk(abs_Cmean, k_level);
CSP_idx = sort(CSP_idx);
end