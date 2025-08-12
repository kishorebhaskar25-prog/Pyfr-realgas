<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%inherit file='base'/>
<%include file='pyfr.solvers.aceuler.kernels.flux'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>

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

    // Compute the flux
    fpdtype_t ftemp[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'R', 'a', 'b', 'cv', 'T', 'pinf')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'{smats}[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
