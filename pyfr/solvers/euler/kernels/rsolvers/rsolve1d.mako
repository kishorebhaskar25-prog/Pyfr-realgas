<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf, R=${R}, a=${a}, b=${b}, cv=${cv}'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];
    fpdtype_t Rgas = R, ag = a, bg = b, cvg = cv;

    utl[0] = ul[0];
    utr[0] = ur[0];
    utl[${nvars - 1}] = ul[${nvars - 1}];
    utr[${nvars - 1}] = ur[${nvars - 1}];

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', off=1)};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', off=1)};

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ntf', 'Rgas', 'ag', 'bg', 'cvg')};

    nf[0] = ntf[0];
    nf[${nvars - 1}] = ntf[${nvars - 1}];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', off=1)};
</%pyfr:macro>
