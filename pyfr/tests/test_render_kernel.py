import importlib.machinery as imm
import importlib.util as imu
import sys
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

# Load makoutil without importing the full pyfr.backends package to avoid
# mpi4py dependency issues during testing.
makoutil_path = Path(__file__).resolve().parents[1] / 'backends' / 'base' / 'makoutil.py'
loader = imm.SourceFileLoader('makoutil_stub', str(makoutil_path))
spec = imu.spec_from_loader(loader.name, loader)
makoutil = imu.module_from_spec(spec)
loader.exec_module(makoutil)
sys.modules['makoutil_stub'] = makoutil


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
${x} ${x}
'''
    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'x': 1 / 3}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert src.count('0.3333333333333333') == 2


def test_render_kernel_helper_macro():
    tplsrc = """
<%namespace module='makoutil_stub' name='pyfr'/>
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${pyfr.dot('a[{i}]', i=3)}
"""
    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    src, *_ = provider._render_kernel('foo', 'foo', {}, {})

    assert 'a[0]' in src and 'a[1]' in src and 'a[2]' in src


def test_macro_numeric_expansion():
    tplsrc = """
<%namespace module='makoutil_stub' name='pyfr'/>
<%
_kernel_argspecs['foo'] = (0, [], [])
import numpy as np
%>
<%pyfr:macro name='coeff' params='i, x'>
idx = i;
val = x;
</%pyfr:macro>
${pyfr.expand('coeff', np.int32(2), np.float64(1/3))}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    src, *_ = provider._render_kernel('foo', 'foo', {}, {})

    assert 'idx = 2;' in src
    assert 'val = 0.3333333333333333;' in src
