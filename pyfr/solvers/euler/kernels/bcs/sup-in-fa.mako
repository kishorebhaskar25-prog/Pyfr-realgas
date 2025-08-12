<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, R, a, b, cv'
             externs='ploc, t'>
    fpdtype_t Rgas = R, ag = a, bg = b, cvg = cv;
    fpdtype_t rho = ${c['rho']};
    fpdtype_t v = 1.0/rho;
    fpdtype_t p = ${c['p']};
    fpdtype_t T = (p + ag/(v*v))*(v - bg)/Rgas;

    ur[0] = rho;
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = rho*(${c[v]});
% endfor
    fpdtype_t e_int = cvg*T - ag*rho;
    ur[${nvars - 1}] = rho*e_int
                     + 0.5*(1.0/rho)*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>
