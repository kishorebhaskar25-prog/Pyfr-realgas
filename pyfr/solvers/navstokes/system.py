from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.navstokes.elements import NavierStokesElements
from pyfr.solvers.navstokes.inters import (NavierStokesBaseBCInters,
                                           NavierStokesIntInters,
                                           NavierStokesMPIInters)
from pyfr.thermo import real_gas as rg


class NavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'navier-stokes'

    elementscls = NavierStokesElements
    intinterscls = NavierStokesIntInters
    mpiinterscls = NavierStokesMPIInters
    bbcinterscls = NavierStokesBaseBCInters

    def __init__(self, backend, mesh, initsoln, nregs, cfg):
        # Ensure real-gas constants are available to all kernels
        rg_consts = [('R', rg.R), ('a', rg.A), ('b', rg.B), ('cv', rg.CV)]
        for k, v in rg_consts:
            if not cfg.hasopt('constants', k):
                cfg.set('constants', k, v)

        super().__init__(backend, mesh, initsoln, nregs, cfg)
