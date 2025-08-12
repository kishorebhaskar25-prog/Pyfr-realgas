<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t rho = ${c['rho']};
    ur[0] = rho;
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = rho*(${c[v]});
% endfor
    fpdtype_t T = (${c['p']} + ${c['a']}*rho*rho)*(1.0 - ${c['b']}*rho)
                /(${c['R']}*rho);
    fpdtype_t E = rho*(${c['cv']}*T - ${c['a']}*rho)
                + 0.5*(1.0/rho)*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
    ur[${nvars - 1}] = E;
</%pyfr:macro>
