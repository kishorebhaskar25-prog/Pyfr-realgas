from collections.abc import Iterable
import itertools as it
import numbers
import re

from mako.runtime import supports_caller, capture
import numpy as np

import pyfr.nputil as nputil
import pyfr.util as util


def ndrange(context, *args):
    for arg in args:
        if isinstance(arg, Iterable):
            if not all(isinstance(a, numbers.Integral) for a in arg):
                raise TypeError(
                    'ndrange expects iterable of integers, got '
                    f'{arg}'
                )
        elif not isinstance(arg, numbers.Integral):
            raise TypeError(
                'ndrange expects integer arguments, got '
                f'{type(arg).__name__}'
            )

    return util.ndrange(*args)


def ilog2range(context, x):
    if not isinstance(x, numbers.Integral):
        raise TypeError(
            'ilog2range expects an integer, got '
            f'{type(x).__name__}'
        )

    return [2**i for i in range(x.bit_length() - 2, -1, -1)]


def npdtype_to_ctype(context, dtype):
    return nputil.npdtype_to_ctype(dtype)


def dot(context, a_, b_=None, /, **kwargs):
    ix, nd = util.first(kwargs.items())
    ab = '({})*({})'.format(a_, b_ or a_)

    # Allow for flexible range arguments
    nd = nd if isinstance(nd, Iterable) else [nd]
    if not all(isinstance(n, numbers.Integral) for n in nd):
        raise TypeError(
            'dot macro expects integer range bounds, got '
            f'{nd}'
        )

    return '(' + ' + '.join(ab.format(**{ix: i}) for i in range(*nd)) + ')'


def array(context, expr_, vals_={}, /, **kwargs):
    ix = util.first(kwargs)
    ni = kwargs.pop(ix)
    items = []

    # Allow for flexible range arguments
    if isinstance(ni, Iterable):
        if not all(isinstance(n, numbers.Integral) for n in ni):
            raise TypeError(
                'array macro expects integer range values for '
                f'{ix}, got {ni}'
            )
        rng = ni
    elif isinstance(ni, numbers.Integral):
        rng = [ni]
    else:
        raise TypeError(
            'array macro expects integer or iterable range for '
            f'{ix}, got {type(ni).__name__}'
        )

    for i in range(*rng):
        if kwargs:
            items.append(array(context, expr_, vals_ | {ix: i}, **kwargs))
        else:
            items.append(expr_.format_map(vals_ | {ix: i}))

    return '{ ' + ', '.join(items) + ' }'


def polyfit(context, f, a, b, n, var, nqpts=500):
    for val, name in [(a, 'a'), (b, 'b')]:
        if not isinstance(val, numbers.Real):
            raise TypeError(
                f'polyfit: {name} must be a real number, got '
                f'{type(val).__name__}'
            )
    for val, name in [(n, 'n'), (nqpts, 'nqpts')]:
        if not isinstance(val, numbers.Integral):
            raise TypeError(
                f'polyfit: {name} must be an integer, got '
                f'{type(val).__name__}'
            )

    x = np.linspace(a, b, nqpts)
    y = f(x)

    coeffs = np.polynomial.polynomial.polyfit(x, y, n)
    pfexpr = f' + {var}*('.join(str(c) for c in coeffs) + ')'*n

    return f'({pfexpr})'


def _strip_parens(s):
    out, depth = [], 0

    for c in s:
        depth += (c in '{(') - (c in ')}')

        if depth == 0 and c not in ')}':
            out.append(c)

    return ''.join(out)


def _locals(body):
    # First, strip away any comments
    body = re.sub(r'//.*?\n', '', body)

    # Next, find all variable declaration statements
    decls = re.findall(r'(?:[A-Za-z_]\w*)\s+([A-Za-z_]\w*[^;]*?);', body)

    # Strip anything inside () or {}
    decls = [_strip_parens(d) for d in decls]

    # A statement can define multiple variables, so split by ','
    decls = it.chain.from_iterable(d.split(',') for d in decls)

    # Extract the variable names
    lvars = [re.match(r'\s*(\w+)', v)[1] for v in decls]

    # Prune invalid names
    return [lv for lv in lvars if lv != 'if']


@supports_caller
def macro(context, name, params, externs=''):
    # Macros may be formed from numeric template arguments; ensure the
    # name is a string prior to any further processing.
    name = str(name)

    # Check we have not already been defined
    if name in context['_macros']:
        raise RuntimeError(f'Attempt to redefine macro "{name}"')

    # Split up the parameter and external variable list
    params = [p.strip() for p in params.split(',')]
    externs = [e.strip() for e in externs.split(',')] if externs else []

    # Ensure no invalid characters in params/extern variables
    for p in it.chain(params, externs):
        if not re.match(r'[A-Za-z_]\w*$', p):
            raise ValueError(f'Invalid param "{p}" in macro "{name}"')

    # Capture the function body
    body = capture(context, context['caller'].body)

    # Identify any local variable declarations
    lvars = _locals(body)

    # Suffix these variables by a '_'
    if lvars:
        body = re.sub(r'\b({0})\b'.format('|'.join(lvars)), r'\1_', body)

    # Save
    context['_macros'][name] = (params, externs, body)

    return ''


def expand(context, name, /, *args, **kwargs):
    # Expand may be called with numeric arguments originating from numpy
    # types.  Coerce these into the corresponding Python primitives so
    # that any downstream formatting, identifier construction or further
    # manipulation behaves as expected.
    def _coerce(v):
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        elif isinstance(v, (int, np.integer)):
            return int(v)
        elif isinstance(v, (float, np.floating)):
            return float(v)
        elif isinstance(v, (str, bytes, np.str_, np.bytes_)):
            try:
                float(v)
            except (TypeError, ValueError):
                return str(v)
            else:
                raise TypeError('numeric value provided as string/bytes')
        elif isinstance(v, list):
            return [_coerce(i) for i in v]
        elif isinstance(v, tuple):
            return tuple(_coerce(i) for i in v)
        elif isinstance(v, dict):
            return {k: _coerce(val) for k, val in v.items()}
        else:
            return v

    # Get the macro parameter list and the body
    name = str(name)
    mparams, mexterns, body = context['_macros'][name]

    # Ensure an appropriate number of arguments have been passed
    if len(mparams) != len(args) + len(kwargs):
        raise ValueError(f'Inconsistent macro parameter list in {name}')

    # Parse the parameter list
    params = dict(zip(mparams, map(_coerce, args)))
    for k, v in kwargs.items():
        if k in params:
            raise ValueError(f'Duplicate macro parameter {k} in {name}')

        params[k] = _coerce(v)

    # Ensure all parameters have been passed
    if sorted(mparams) != sorted(params):
        raise ValueError(f'Inconsistent macro parameter list in {name}')

    # Ensure all (used) external parameters have been passed to the kernel
    for extrn in mexterns:
        if (extrn not in context['_extrns'] and
            re.search(rf'\b{extrn}\b', body)):
            raise ValueError(f'Missing external {extrn} in {name}')

    # Rename local parameters
    for lname, subst in params.items():
        body = re.sub(rf'\b{lname}\b', str(subst), body)

    return f'{{\n{body}\n}}'


@supports_caller
def kernel(context, name, ndim, **kwargs):
    # Kernel names may be formed dynamically from numeric template
    # arguments.  Ensure the resulting name is a string so that it can be
    # used as a dictionary key and passed to the generator class.
    name = str(name)

    extrns = context['_extrns']

    # Validate the argument list
    if any(arg in extrns for arg in kwargs):
        raise ValueError('Duplicate argument in {0}: {1} {2}'
                         .format(name, list(kwargs), list(extrns)))

    # Merge local and external arguments
    kwargs = dict(kwargs, **extrns)

    # Capture the kernel body
    body = capture(context, context['caller'].body)

    # Get the generator class and data types
    kerngen = context['_kernel_generator']
    fpdtype, ixdtype = context['fpdtype'], context['ixdtype']

    # Instantiate
    kern = kerngen(name, int(ndim), kwargs, body, fpdtype, ixdtype)

    # Save the argument/type list for later use
    context['_kernel_argspecs'][name] = kern.argspec()

    # Render and return the complete kernel
    return kern.render()


def alias(context, name, func):
    # Allow for macros to be aliased using numeric template arguments.
    name = str(name)
    func = str(func)

    context['_macros'][name] = context['_macros'][func]
    return ''
