<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.aceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.acnavstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.acnavstokes.kernels.flux'/>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur, Rgas, ag, bg, cvg, cpg, T0, p0'>
    fpdtype_t nor = ${' + '.join(f'nl[{i}]*ul[{i + 1}]' for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, nl, magnl, Rgas, ag, bg, cvg, cpg, T0, p0'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur', 'Rgas', 'ag', 'bg', 'cvg', 'cpg', 'T0', 'p0')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*ficomm[${i}];
% endfor
</%pyfr:macro>
