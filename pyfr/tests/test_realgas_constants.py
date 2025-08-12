from importlib import machinery as imm, util as imu
from pathlib import Path

from mako.template import Template

from pyfr.inifile import Inifile
from pyfr.thermo import real_gas as rg


# Import BasePointwiseKernelProvider without initialising a backend
kernel_path = Path(__file__).resolve().parents[1] / 'backends' / 'base' / 'kernels.py'
loader = imm.SourceFileLoader('kernels', str(kernel_path))
spec = imu.spec_from_loader(loader.name, loader)
kernels = imu.module_from_spec(spec)
loader.exec_module(kernels)
BasePointwiseKernelProvider = kernels.BasePointwiseKernelProvider


class DummyLookup:
    def __init__(self, src):
        self.src = src

    def get_template(self, mod):
        return Template(self.src)


class DummyBackend:
    def __init__(self, src):
        self.lookup = DummyLookup(src)


class DummyProvider(BasePointwiseKernelProvider):
    kernel_generator_cls = object


def test_constants_arrive_as_floats_in_kernel():
    cfg = Inifile()

    # Populate the standard real-gas constants
    for k, v in [('R', rg.R), ('a', rg.A), ('b', rg.B), ('cv', rg.CV)]:
        if not cfg.hasopt('constants', k):
            cfg.set('constants', k, v)

    consts = cfg.items_as('constants', float)
    assert all(isinstance(v, float) for v in consts.values())

    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${c['R'] + c['a'] + c['b'] + c['cv']}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    src, *_ = provider._render_kernel('foo', 'foo', {}, {'c': consts})

    # The rendered source should contain the numeric sum of the constants
    expected = sum(consts.values())
    assert str(expected) in src
