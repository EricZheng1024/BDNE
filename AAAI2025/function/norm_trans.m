% function objs = norm_trans(objs,zmin,zmax,alpha,gain)
%     objs_n = (objs-zmin)*gain ./ (zmax-zmin);
%     objs = (1-alpha)*objs_n + alpha*mean(objs_n,1);  % 列向量
% end
function objs = norm_trans(objs,zmin,zmax,alpha_gain,minus_alpha_gain,gain)
    objs_n = (objs-zmin)*gain ./ (zmax-zmin);
    objs = minus_alpha_gain*objs_n + alpha_gain*mean(objs_n,1);  % 列向量
end


% function objs_n = normalize_2(objs,zmin,zmax)
%     objs_n = (objs-zmin) ./ (zmax-zmin);
% end


% function objs = trans_cone(objs,alpha)
%     objs = (1-alpha)*objs + alpha*mean(objs,1);  % 列向量
% end