<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t inv = 1.0/ul[0];
    fpdtype_t e_i = ul[${nvars - 1}] - 0.5*inv*
                    ${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t T_i = (e_i*inv + ${c['a']}*ul[0])/${c['cv']};
    fpdtype_t pl = ${c['R']}*T_i*ul[0]/(1.0 - ${c['b']}*ul[0])
                 - ${c['a']}*ul[0]*ul[0];
    fpdtype_t Rdcp = ${c['R']}/(${c['cv']} + ${c['R']});
    fpdtype_t udotu = ${2.0*c['cpTt']}*(1.0
                    - pow(${c['pt']}, -Rdcp)*pow(pl, Rdcp));
    udotu = fmax(0, udotu);

    fpdtype_t rho = (${c['cv']} + ${c['R']})*pl/(${c['R']}*(${c['cpTt']} - 0.5*udotu));
    ur[0] = rho;
% for i, v in enumerate(c['vc']):
    ur[${i + 1}] = ${v}*rho*sqrt(udotu);
% endfor
    fpdtype_t E = (${c['cv']}/(${c['cv']} + ${c['R']}))*rho*(${c['cpTt']} - 0.5*udotu)
                - ${c['a']}*rho*rho + 0.5*rho*udotu;
    ur[${nvars - 1}] = E;
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
