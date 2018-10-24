% GEO_FUNCTIONS
%   Stores in geo_exp, geo_log, geo_norm, geo_map, geo_tp and geo_ps the
%   correct function as a handle.

% This file comes from the project "C1 bezier paths on surfaces"
% by Gousenbourger et al. 
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt


switch manifold
    case 'euclidean'
        geo_exp  = @exp_eucl;
        geo_log  = @log_eucl;
        geo_norm = @norm_eucl;
        geo_map  = @(x,s,t,m) exp_eucl(x,log_eucl(x,s),t);
        geo_tp   = @tp_eucl;
        geo_ps   = @ps_eucl;
        geo_dist = @dist_eucl_geo;
    case 'sphere'
        geo_exp  = @exp_sphere;
        geo_log  = @log_sphere;
        geo_norm = @norm_sphere;
        geo_map  = @(x,s,t,m) exp_sphere(x,log_sphere(x,s),t);
        geo_tp   = @tp_sphere;
        geo_ps   = @ps_sphere;
        geo_dist = @dist_sphere_geo;
    case 'psd_quotient'
        geo_exp  = @exp_psd_quotient;
        geo_log  = @log_psd_quotient;
        geo_norm = @norm_psd_quotient;
        geo_map  = @(x,s,t,m) exp_psd_quotient(x,log_psd_quotient(x,s),t);
        geo_tp   = @tp_psd_quotient;
        geo_ps   = @ps_psd_quotient;
        geo_dist = @dist_psd_quotient_geo;
    case 'SO3'
        geo_exp  = @exp_so3;
        geo_log  = @log_so3;
        geo_norm = @norm_so3;
        geo_map  = @(x,s,t,m) exp_so3(x,log_so3(x,s),t);
        geo_tp   = @tp_so3;
        geo_ps   = @ps_so3;
        geo_dist = @dist_so3_geo;
    otherwise error('Wrong manifold')
end
