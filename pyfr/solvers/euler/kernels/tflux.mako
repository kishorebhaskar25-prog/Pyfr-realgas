<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<% smats = 'smats_l' if 'linear' in ktype else 'smats' %>

<%pyfr:kernel name='tflux' ndim='2'
              u='in fpdtype_t[${fmt(nvars)}]'
              f='out fpdtype_t[${fmt(ndims)}][${fmt(nvars)}]'
              smats='in fpdtype_t[${fmt(ndims)}][${fmt(ndims)}]'
              verts='in broadcast-col fpdtype_t[${fmt(nverts)}][${fmt(ndims)}]'
              upts='in broadcast-row fpdtype_t[${fmt(ndims)}]'>
% if 'linear' in ktype:
    // Compute the S matrices
    fpdtype_t ${smats}[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', smats, 'djac')};
% endif

    // Compute the flux using the real-gas equation of state
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t p, c, v[${ndims}];
    fpdtype_t Rgas = ${c['R']}, ag = ${c['a']}, bg = ${c['b']}, cvg = ${c['cv']};
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'p', 'c', 'v', 'Rgas', 'ag', 'bg', 'cvg')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'{smats}[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
