from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.euler.elements import EulerElements
from pyfr.solvers.euler.inters import (EulerIntInters, EulerMPIInters,
                                       EulerBaseBCInters)
from pyfr.thermo import real_gas as rg


class EulerSystem(BaseAdvectionSystem):
    name = 'euler'

    elementscls = EulerElements
    intinterscls = EulerIntInters
    mpiinterscls = EulerMPIInters
    bbcinterscls = EulerBaseBCInters

    def __init__(self, backend, mesh, initsoln, nregs, cfg):
        # Ensure real-gas constants are available to all kernels
        for k, v in [('R', rg.R), ('a', rg.A), ('b', rg.B), ('cv', rg.CV)]:
            if not cfg.hasopt('constants', k):
                cfg.set('constants', k, str(v))

        super().__init__(backend, mesh, initsoln, nregs, cfg)
