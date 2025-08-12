<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%inherit file='base'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${fmt(nvars)}]'
              gradul='in view fpdtype_t[${fmt(ndims)}][${fmt(nvars)}]'
              artviscl='in view fpdtype_t'
              nl='in fpdtype_t[${fmt(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv};
    ${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'artviscl', 'norm_nl', 'mag_nl', 'Rgas', 'ag', 'bg', 'cvg')};
</%pyfr:kernel>
