<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% gamma = 1.0 + c['R']/c['cv'] %>
<% gmo = gamma - 1.0 %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t rhoe = ${c['rho']};
    fpdtype_t p0 = ${c['p']};
    fpdtype_t one_m_br = 1.0 - ${c['b']}*rhoe;
    fpdtype_t Te = (p0 + ${c['a']}*rhoe*rhoe)*one_m_br/(${c['R']}*rhoe);
    fpdtype_t cs = sqrt(${c['R']}*Te*(1.0 + ${c['R']}/${c['cv']})/(one_m_br*one_m_br)
                      - 2*${c['a']}*rhoe);
    fpdtype_t s = (p0 + ${c['a']}*rhoe*rhoe)*pow(rhoe, -${gamma})
                  *pow(one_m_br, ${gamma});
    fpdtype_t ratio = cs*${2.0/gmo};

    fpdtype_t inv = 1.0/ul[0];
    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};
    fpdtype_t V_i = inv*(${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                      for i in range(ndims))});
    fpdtype_t e_i = ul[${nvars - 1}]
                  - 0.5*inv*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t T_i = (e_i*inv + ${c['a']}*ul[0])/${c['cv']};
    fpdtype_t p_i = ${c['R']}*T_i*ul[0]/(1.0 - ${c['b']}*ul[0])
                  - ${c['a']}*ul[0]*ul[0];
    fpdtype_t c_i = sqrt(${c['R']}*T_i*(1.0 + ${c['R']}/${c['cv']})/
                          ((1.0 - ${c['b']}*ul[0])*(1.0 - ${c['b']}*ul[0]))
                          - 2*${c['a']}*ul[0]);
    fpdtype_t R_e = (fabs(V_e) >= cs && V_i >= 0)
                  ? V_i - c_i*${2.0/gmo}
                  : V_e - ratio;
    fpdtype_t R_i = (fabs(V_e) >= cs && V_i < 0)
                  ? V_e + ratio
                  : V_i + c_i*${2.0/gmo};
    fpdtype_t V_b = 0.5*(R_e + R_i);
    fpdtype_t c_b = ${0.25*gmo}*(R_i - R_e);
    fpdtype_t rho_b = (V_i < 0)
                    ? pow((1.0/(${gamma}*s))*c_b*c_b, ${1.0/gmo})
                    : ul[0]*pow(ul[0]*c_b*c_b/(${gamma}*p_i), ${1.0/gmo});
    fpdtype_t p_b = ${c['R']}*Te*pow(rho_b/rhoe, ${gamma})/(1.0 - ${c['b']}*rho_b)
                   - ${c['a']}*rho_b*rho_b;

    ur[0] = rho_b;
% for i in range(ndims):
    ur[${i + 1}] = (V_i >= 0)
                 ? rho_b*(ul[${i + 1}]*inv + (V_b - V_i)*nl[${i}])
                 : rho_b*(${c['uvw'[i]]} + (V_b - V_e)*nl[${i}]);
% endfor
    ur[${nvars - 1}] = p_b*${1.0/gmo}
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>
