import importlib.machinery as imm
import importlib.util as imu
from pathlib import Path

import numpy as np
from mako.template import Template

# Import BasePointwiseKernelProvider without triggering backend initialisation
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


def test_render_kernel_numeric_coercion():
    tplsrc = '''
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${x} ${y}
'''
    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'x': 1/3, 'y': np.float64(1/7)}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert '0.33333333333333331' in src
    assert '0.14285714285714285' in src
