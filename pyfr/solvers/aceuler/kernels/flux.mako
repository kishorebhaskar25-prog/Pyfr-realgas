<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%!
from pyfr.thermo import real_gas as rg
%>

<%pyfr:macro name='inviscid_flux' params='s, f, Rgas, ag, bg, cvg, T0, p0'>
    // Base-state density from the van der Waals equation of state
    const fpdtype_t rho0 = ${rg.rho_from_pT(p0, T0, Rgas, ag, bg)};
    const fpdtype_t ombr0 = 1.0 - bg*rho0;
    const fpdtype_t p0v = (rho0*Rgas*T0)/ombr0 - ag*rho0*rho0;
    const fpdtype_t c0 = sqrt((Rgas*T0/(ombr0*ombr0))*(1.0 + Rgas/cvg)
                              - 2.0*ag*rho0);
    (void)p0v; (void)c0;

    // Velocity components
    fpdtype_t v[] = ${pyfr.array('s[{i}]', i=(1, ndims + 1))};

    // Pressure perturbation
    fpdtype_t p = s[0];

    // Mass flux (rho0 * v)
% for i in range(ndims):
    f[${i}][0] = rho0*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = v[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
