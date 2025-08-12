<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, R, a, b, cv' externs='ploc, t'>
    fpdtype_t Rgas = R, ag = a, bg = b, cvg = cv;
% for i in range(nvars - 1):
    ur[${i}] = ul[${i}];
% endfor
    fpdtype_t rho = ur[0], invrho = 1.0/rho;
    fpdtype_t p = ${c['p']};
    fpdtype_t T = (p + ag/(invrho*invrho))*(invrho - bg)/Rgas;
    fpdtype_t e_int = cvg*T - ag*rho;
    ur[${nvars - 1}] = rho*e_int
                     + 0.5*invrho*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
