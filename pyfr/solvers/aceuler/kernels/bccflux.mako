<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.aceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.aceuler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${fmt(nvars)}]'
              nl='in fpdtype_t[${fmt(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv},
              T0 = ${T}, p0 = ${pinf};

    // Compute the RHS
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'norm_nl', 'ur',
                  'Rgas', 'ag', 'bg', 'cvg', 'T0', 'p0')};

    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'fn', 'Rgas', 'ag', 'bg', 'cvg', 'T0', 'p0')};

    // Scale and write out the common normal fluxes
% for i in range(nvars):
    ul[${i}] = mag_nl*fn[${i}];
% endfor
</%pyfr:kernel>
