<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.acnavstokes.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bcconu' ndim='1'
              ulin='in view fpdtype_t[${fmt(nvars)}]'
              ulout='out view fpdtype_t[${fmt(nvars)}]'
              nlin='in fpdtype_t[${fmt(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nlin[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nlin[{i}]', i=ndims)};

    fpdtype_t Rgas = ${R}, ag = ${a}, bg = ${b}, cvg = ${cv},
              T0 = ${T}, p0 = ${pinf};

    ${pyfr.expand('bc_ldg_state', 'ulin', 'norm_nl', 'ulout',
                  'Rgas', 'ag', 'bg', 'cvg', 'T0', 'p0')};
</%pyfr:kernel>
