from pyfr.solvers.aceuler.elements import BaseACFluidElements
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.thermo import real_gas as rg


class ACNavierStokesElements(BaseACFluidElements,
                             BaseAdvectionDiffusionElements):
    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        return grad_cons

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        kprefix = 'pyfr.solvers.acnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.tflux')

        # Real-gas constants
        c = self.cfg.items_as('constants', float)
        c.setdefault('R', rg.R)
        c.setdefault('a', rg.A)
        c.setdefault('b', rg.B)
        c.setdefault('cv', rg.CV)
        c.setdefault('cp', c['cv'] + c['R'])
        c.setdefault('T', 0.0)
        c.setdefault('pinf', 0.0)

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': c,
            'jac_exprs': self.basis.jac_exprs,
            'R': c['R'], 'a': c['a'], 'b': c['b'],
            'cv': c['cv'], 'cp': c['cp'], 'T': c['T'],
            'pinf': c['pinf']
        }

        # Helpers
        tdisf = []
        crv, lin = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        # Gradient + flux kernel fusion
        if self.grad_fusion:
            if crv in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved-fused'},
                    dims=[self.nupts, r[crv]],  u=s(self.scal_upts[uin], crv),
                    f=s(self._vect_upts, crv), gradu=s(self._grad_upts, crv),
                    rcpdjac=self.rcpdjac_at('upts', 'curved'),
                    smats=self.curved_smat_at('upts')
                ))
            if lin in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear-fused'},
                    dims=[self.nupts, r[lin]], u=s(self.scal_upts[uin], lin),
                    f=s(self._vect_upts, lin), gradu=s(self._grad_upts, lin),
                    verts=self.ploc_at('linspts', lin), upts=self.upts
                ))

            def tdisf_k(uin):
                return self._make_sliced_kernel(k(uin) for k in tdisf)

            self.kernels['tdisf_fused'] = tdisf_k
        # No gradient + flux kernel fusion, with flux-AA
        elif 'flux' in self.antialias:
            if crv in r:
                tdisf.append(lambda: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved'},
                    dims=[self.nqpts, r[crv]], u=s(self._scal_qpts, crv),
                    f=s(self._vect_qpts, crv), smats=self.curved_smat_at('qpts')
                ))
            if lin in r:
                tdisf.append(lambda: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear'},
                    dims=[self.nqpts, r[lin]], u=s(self._scal_qpts, lin),
                    f=s(self._vect_qpts, lin), verts=self.ploc_at('linspts', lin),
                    upts=self.qpts
                ))

            def tdisf_k():
                return self._make_sliced_kernel(k() for k in tdisf)

            self.kernels['tdisf'] = tdisf_k
        # No gradient + flux kernel fusion, no flux-AA
        else:
            if crv in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved'},
                    dims=[self.nupts, r[crv]], u=s(self.scal_upts[uin], crv),
                    f=s(self._vect_upts, crv), smats=self.curved_smat_at('upts')
                ))
            if lin in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear'},
                    dims=[self.nupts, r[lin]], u=s(self.scal_upts[uin], lin),
                    f=s(self._vect_upts, lin), verts=self.ploc_at('linspts', lin),
                    upts=self.upts
                ))

            def tdisf_k(uin):
                return self._make_sliced_kernel(k(uin) for k in tdisf)

            self.kernels['tdisf'] = tdisf_k
