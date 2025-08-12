<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, Rgas, ag, bg, cvg, cpg, T0, p0' externs='ploc, t'>
    ur[0] = ${c['p']};

% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = ul[${i + 1}];
% endfor
</%pyfr:macro>
