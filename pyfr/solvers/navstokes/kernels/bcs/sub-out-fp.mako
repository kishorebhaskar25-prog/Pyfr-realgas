<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nvars - 1):
    ur[${i}] = ul[${i}];
% endfor
    fpdtype_t rho = ul[0];
    fpdtype_t T = (${c['p']} + ${c['a']}*rho*rho)*(1.0 - ${c['b']}*rho)
                /(${c['R']}*rho);
    fpdtype_t E = rho*(${c['cv']}*T - ${c['a']}*rho)
                + 0.5*(1.0/rho)*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
    ur[${nvars - 1}] = E;
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
