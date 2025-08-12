import importlib.machinery as imm
import importlib.util as imu
from pathlib import Path

import numpy as np
from mako.template import Template


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


def test_float_coercion():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
name = f'coef_{x}'
%>
${name}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'x': np.float64(1 / 3)}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert 'coef_0.3333333333333333' in src


def test_integer_range_loop():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
% for i in range(n):
${i}
% endfor
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'n': 4}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert '0\n1\n2\n3' in src


def test_numeric_param_multiplication():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${a*b}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'a': np.int32(2), 'b': np.float64(1 / 3)}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert '0.6666666666666666' in src


def test_string_tplarg_accepted():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${x}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, {'x': '1'})

    assert '1' in src


def test_bytes_tplarg_accepted():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
${x}
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, {'x': b'1'})

    assert "b'1'" in src


def test_string_list_tplarg_accepted():
    tplsrc = """
<%
_kernel_argspecs['foo'] = (0, [], [])
%>
% for e in exprs:
${e}
% endfor
"""

    backend = DummyBackend(tplsrc)
    provider = DummyProvider(backend)

    tplargs = {'exprs': ['x + 1', 'y - 2']}
    src, ndim, argn, argt = provider._render_kernel('foo', 'foo', {}, tplargs)

    assert 'x + 1' in src and 'y - 2' in src

