<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities, pressures and sounds
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, cl, cr;

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv};

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'cl', 'vl', 'Rgas', 'ag', 'bg', 'cvg')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'cr', 'vr', 'Rgas', 'ag', 'bg', 'cvg')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    fpdtype_t a = 0.5*(cl + cr) + 0.25*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
