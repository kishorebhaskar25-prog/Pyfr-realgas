<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%!
from pyfr.thermo import real_gas as rg

def pressure(rho, T, R, a, b):
    return f'(({R}*{T})/(1.0/{rho} - {b}) - {a}/((1.0/{rho})*(1.0/{rho})))'

def energy(rho, T, cv, a):
    return f'({cv}*{T} - {a}*{rho})'

def T_from_rho_e(rho, e, cv, a):
    return f'(({e} + {a}*{rho})/{cv})'
%>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, R, a, b, cv'
             externs='ploc, t, ic, im'>
    fpdtype_t Rgas = R, ag = a, bg = b, cvg = cv;
    fpdtype_t gmo = Rgas/cvg;

    fpdtype_t p = t*im + ic;
    fpdtype_t rho_e = ${c['rho']};
    fpdtype_t v_e = 1.0/rho_e;
    fpdtype_t T_e = (p + ag/(v_e*v_e))*(v_e - bg)/Rgas;
    fpdtype_t ombr_e = 1.0 - bg*rho_e;
    fpdtype_t cs = sqrt((Rgas*T_e/(ombr_e*ombr_e))*(1.0 + gmo)
                        - 2.0*ag*rho_e);
    fpdtype_t s = p*pow(rho_e, -(1.0 + gmo));
    fpdtype_t ratio = cs*(2.0/gmo);

    fpdtype_t inv = 1.0/ul[0];
    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};
    fpdtype_t V_i = inv*(${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                      for i in range(ndims))});
    fpdtype_t E_i = ul[${nvars - 1}];
    fpdtype_t e_i = E_i - 0.5*inv*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t T_i = ${T_from_rho_e('ul[0]', 'e_i', 'cvg', 'ag')};
    fpdtype_t p_i = ${pressure('ul[0]', 'T_i', 'Rgas', 'ag', 'bg')};
    fpdtype_t ombr_i = 1.0 - bg*ul[0];
    fpdtype_t c_i = sqrt((Rgas*T_i/(ombr_i*ombr_i))*(1.0 + gmo)
                         - 2.0*ag*ul[0]);
    fpdtype_t R_e = (fabs(V_e) >= cs && V_i >= 0)
                  ? V_i - c_i*(2.0/gmo)
                  : V_e - ratio;
    fpdtype_t R_i = (fabs(V_e) >= cs && V_i < 0)
                  ? V_e + ratio
                  : V_i + c_i*(2.0/gmo);
    fpdtype_t V_b = 0.5*(R_e + R_i);
    fpdtype_t c_b = 0.25*gmo*(R_i - R_e);
    fpdtype_t rho_b = (V_i < 0)
                    ? pow((1.0/((1.0 + gmo)*s))*c_b*c_b, cvg/Rgas)
                    : ul[0]*pow(ul[0]*c_b*c_b/((1.0 + gmo)*p_i), cvg/Rgas);
    fpdtype_t invrho_b = 1.0/rho_b;
    fpdtype_t ombr_b = 1.0 - bg*rho_b;
    fpdtype_t T_b = (c_b*c_b + 2.0*ag*rho_b)*(ombr_b*ombr_b)
                    /(Rgas*(1.0 + gmo));
    fpdtype_t p_b = ${pressure('rho_b', 'T_b', 'Rgas', 'ag', 'bg')};

    ur[0] = rho_b;
% for i in range(ndims):
    ur[${i + 1}] = (V_i >= 0)
                 ? rho_b*(ul[${i + 1}]*inv + (V_b - V_i)*nl[${i}])
                 : rho_b*(${c['uvw'[i]]} + (V_b - V_e)*nl[${i}]);
% endfor
    fpdtype_t e_b = ${energy('rho_b', 'T_b', 'cvg', 'ag')};
    ur[${nvars - 1}] = rho_b*e_b
                     + 0.5*invrho_b*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>
