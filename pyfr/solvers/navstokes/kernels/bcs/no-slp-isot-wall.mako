<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, R, a, b, cv' externs='ploc, t'>
    fpdtype_t Rgas = R, ag = a, cvg = cv;
    fpdtype_t cpT = ${c['cpTw']};
    fpdtype_t Tw = cpT/(cvg + Rgas);

    ur[0] = ul[0];
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
% endfor
    fpdtype_t e_int = cvg*Tw - ag*ur[0];
    ur[${nvars - 1}] = ur[0]*e_int
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur, R, a, b, cv' externs='ploc, t'>
    fpdtype_t Rgas = R, ag = a, cvg = cv;
    fpdtype_t cpT = ${c['cpTw']};
    fpdtype_t Tw = cpT/(cvg + Rgas);

    ur[0] = ul[0];
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = ${c[v]}*ul[0];
% endfor
    fpdtype_t e_int = cvg*Tw - ag*ur[0];
    ur[${nvars - 1}] = ur[0]*e_int
                     + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
