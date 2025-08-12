from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.thermo import real_gas as rg


class BaseACFluidElements:
    @staticmethod
    def privars(ndims, cfg):
        return ['p', 'u', 'v'] if ndims == 2 else ['p', 'u', 'v', 'w']

    convars = privars

    @staticmethod
    def dualcoeffs(ndims, cfg):
        return ['u', 'v'] if ndims == 2 else ['u', 'v', 'w']

    @staticmethod
    def visvars(ndims, cfg):
        if ndims == 2:
            return {
                'velocity': ['u', 'v'],
                'pressure': ['p']
            }
        elif ndims == 3:
            return {
                'velocity': ['u', 'v', 'w'],
                'pressure': ['p']
            }

    @staticmethod
    def pri_to_con(pris, cfg):
        return list(pris)

    @staticmethod
    def con_to_pri(convs, cfg):
        return list(convs)

    @staticmethod
    def diff_con_to_pri(cons, diff_cons, cfg):
        return list(diff_cons)

    @staticmethod
    def validate_formulation(controller):
        if controller.formulation != 'dual':
            raise ValueError('System not compatible with time stepping '
                             'formulation.')


class ACEulerElements(BaseACFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.aceuler.kernels.tflux')

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
        slicedk = self._make_sliced_kernel

        if crv in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nupts, r[crv]], u=s(self.scal_upts[uin], crv),
                f=s(self._vect_upts, crv), smats=self.curved_smat_at('upts')
            ))
        elif crv in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nqpts, r[crv]],
                u=s(self._scal_qpts, crv), f=s(self._vect_qpts, crv),
                smats=self.curved_smat_at('qpts')
            ))

        if lin in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nupts, r[lin]], u=s(self.scal_upts[uin], lin),
                f=s(self._vect_upts, lin), verts=self.ploc_at('linspts', lin),
                upts=self.upts
            ))
        elif lin in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nqpts, r[lin]], u=s(self._scal_qpts, lin),
                f=s(self._vect_qpts, lin), verts=self.ploc_at('linspts', lin),
                upts=self.qpts
            ))

        self.kernels['tdisf'] = lambda: slicedk(k() for k in tdisf)
