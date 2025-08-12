<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.acnavstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.acnavstokes.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv},
              T0 = ${T}, p0 = ${pinf};

    ${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'norm_nl', 'mag_nl',
                  'Rgas', 'ag', 'bg', 'cvg', 'T0', 'p0')};
</%pyfr:kernel>
