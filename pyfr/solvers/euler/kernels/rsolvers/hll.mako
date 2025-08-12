<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf, R=${R}, a=${a}, b=${b}, cv=${cv}'>
    // Compute the left and right fluxes + velocities, pressures and sounds
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, al, ar;
    fpdtype_t nf_sub, nf_fl, nf_fr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'al', 'vl', 'R', 'a', 'b', 'cv')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'ar', 'vr', 'R', 'a', 'b', 'cv')};

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Average normal velocity and sound speed
    fpdtype_t nv = (sqrt(ul[0])*nvl + sqrt(ur[0])*nvr)
                 / (sqrt(ul[0]) + sqrt(ur[0]));
    fpdtype_t a = 0.5*(al + ar);

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = min(nv - a, nvl - al);
    fpdtype_t sr = max(nv + a, nvr + ar);
    fpdtype_t rcpsrsl = 1/(sr - sl);

    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join(f'n[{j}]*fl[{j}][{i}]' for j in range(ndims))};
    nf_fr = ${' + '.join(f'n[{j}]*fr[{j}][{i}]' for j in range(ndims))};
    nf_sub = (sr*nf_fl - sl*nf_fr + sl*sr*(ur[${i}] - ul[${i}]))*rcpsrsl;
    nf[${i}] = (0 <= sl) ? nf_fl : (0 >= sr) ? nf_fr : nf_sub;
% endfor
</%pyfr:macro>
