<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.aceuler.kernels.flux'/>

<%!
from pyfr.thermo import real_gas as rg
%>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf, R, a, b, cv, T0, p0'>
    const fpdtype_t rho0 = ${rg.rho_from_pT(p0, T0, R, a, b)};
    const fpdtype_t c0 = sqrt(${rg.sound_speed_sq(rho0, T0, R, a, b, cv)});

    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'R', 'a', 'b', 'T0', 'p0')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'R', 'a', 'b', 'T0', 'p0')};

    fpdtype_t vl[${ndims}] = ${pyfr.array('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t vr[${ndims}] = ${pyfr.array('ur[{i}]', i=(1, ndims + 1))};

    // Normal of the average interface velocity
    fpdtype_t nv = 0.5*${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the wave speed
    fpdtype_t a = fabs(nv) + c0;

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
