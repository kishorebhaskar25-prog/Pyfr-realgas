<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, p, v'>
    fpdtype_t rho = s[0], invrho = 1.0/rho, E = s[${nvars - 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Internal energy and temperature
    fpdtype_t e = E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};
    fpdtype_t T = (e + ${c['a']}*rho)/${c['cv']};
    fpdtype_t rb = 1.0 - ${c['b']}*rho;

    // Compute the pressure from the real-gas equation of state
    p = ${c['R']}*T*rho/rb - ${c['a']}*rho*rho;

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];
    f[${i}][${nvars - 1}] = (E + p)*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='s, f, p, v'>
    fpdtype_t rho = s[0], invrho = 1.0/rho, E = s[${nvars - 1}];

    // Compute the velocities
% for i in range(ndims):
    v[${i}] = invrho*s[${i + 1}];
% endfor

    // Internal energy and temperature
    fpdtype_t e = E - 0.5*invrho*${pyfr.dot('s[{i}]', i=(1, ndims + 1))};
    fpdtype_t T = (e + ${c['a']}*rho)/${c['cv']};
    fpdtype_t rb = 1.0 - ${c['b']}*rho;

    // Compute the pressure from the real-gas equation of state
    p = ${c['R']}*T*rho/rb - ${c['a']}*rho*rho;

    // Density and energy fluxes
    f[0] = s[1];
    f[${nvars - 1}] = (E + p)*v[0];

    // Momentum fluxes
    f[1] = s[1]*v[0] + p;
% for j in range(1, ndims):
    f[${j + 1}] = s[1]*v[${j}];
% endfor
</%pyfr:macro>
