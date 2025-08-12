<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_entropy' params='u, d, p, e, R, a, b, cv'>
    d = u[0];
    fpdtype_t rcpd = 1.0/d;
    fpdtype_t E = u[${nvars - 1}];

    // Internal energy and temperature
    fpdtype_t ei = E - 0.5*rcpd*(${pyfr.dot('u[{i}]', i=(1, ndims + 1))});
    fpdtype_t T = (ei + a*d)/cv;

    // Compute the pressure
    fpdtype_t vvol = rcpd;
    p = (R*T)/(vvol - b) - a/(vvol*vvol);

    // Compute entropy
    fpdtype_t gamma = 1.0 + R/cv;
    e = (d > 0 && p > 0) ? p*pow(rcpd, gamma) : ${fpdtype_max};
</%pyfr:macro>
