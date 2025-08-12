<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Normal velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Sound speeds for left and right states
    fpdtype_t el = ul[${nvars - 1}] - 0.5*(1.0/ul[0])*${pyfr.dot('ul[{i + 1}]', 'ul[{i + 1}]', i=ndims)};
    fpdtype_t Tl = (el + ${c['a']}*ul[0])/${c['cv']};
    fpdtype_t bl = 1.0 - ${c['b']}*ul[0];
    fpdtype_t cl = sqrt(${c['R']}*Tl/(bl*bl) - 2*${c['a']}*ul[0]
                         - (${c['R']}*${c['R']}*Tl)/(${c['cv']}*bl*bl));

    fpdtype_t er = ur[${nvars - 1}] - 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i + 1}]', 'ur[{i + 1}]', i=ndims)};
    fpdtype_t Tr = (er + ${c['a']}*ur[0])/${c['cv']};
    fpdtype_t br = 1.0 - ${c['b']}*ur[0];
    fpdtype_t cr = sqrt(${c['R']}*Tr/(br*br) - 2*${c['a']}*ur[0]
                         - (${c['R']}*${c['R']}*Tr)/(${c['cv']}*br*br));

    fpdtype_t smax = fmax(fabs(nvl) + cl, fabs(nvr) + cr);
    fpdtype_t a = 0.5*smax;

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
