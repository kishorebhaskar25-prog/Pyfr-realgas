<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.bcs.${bctype}'/>
<%include file='pyfr.solvers.euler.kernels.entropy'/>

<%pyfr:kernel name='bccent' ndim='1'
              ul='in view fpdtype_t[${fmt(nvars)}]'
              nl='in fpdtype_t[${fmt(ndims)}]'
              entmin_lhs='out view reduce(min) fpdtype_t'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv};

    // Compute the RHS
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'norm_nl', 'ur', 'Rgas', 'ag', 'bg', 'cvg')};

    // Compute entropy for boundary state
    fpdtype_t p, d, entmin_rhs;
    ${pyfr.expand('compute_entropy', 'ur', 'd', 'p', 'entmin_rhs', 'Rgas', 'ag', 'bg', 'cvg')};

    // Compute face minima (reduce with atomics)
    entmin_lhs = entmin_rhs;
</%pyfr:kernel>
