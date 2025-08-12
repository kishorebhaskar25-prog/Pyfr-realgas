<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur, R, a, b, cv' externs='ploc, t'>
    fpdtype_t Rgas = R, ag = a, bg = b, cvg = cv;
    fpdtype_t cp = cvg + Rgas;
    fpdtype_t Rdcp = Rgas/cp;

    fpdtype_t invrho = 1.0/ul[0];
    fpdtype_t e = ul[${nvars - 1}] - 0.5*invrho*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t T = (e + ag*ul[0])/cvg;
    fpdtype_t pl = (Rgas*T)/(invrho - bg) - ag/(invrho*invrho);

    fpdtype_t cpTt = ${c['cpTt']};
    fpdtype_t pt = ${c['pt']};
    fpdtype_t udotu = 2.0*cpTt*(1.0 - pow(pl/pt, Rdcp));
    udotu = fmax(0, udotu);

    T = (cpTt - 0.5*udotu)/cp;
    fpdtype_t rho = pl/(Rgas*T);
    for (int _i = 0; _i < 5; ++_i) {
        fpdtype_t f = (Rgas*T*rho)/(1.0 - bg*rho) - ag*rho*rho - pl;
        fpdtype_t df = (Rgas*T)/((1.0 - bg*rho)*(1.0 - bg*rho)) - 2.0*ag*rho;
        rho -= f/df;
    }
    ur[0] = rho;
% for i, v in enumerate(c['vc']):
    ur[${i + 1}] = ${v}*rho*sqrt(udotu);
% endfor
    fpdtype_t e_int = cvg*T - ag*rho;
    ur[${nvars - 1}] = rho*e_int + 0.5*rho*udotu;
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
