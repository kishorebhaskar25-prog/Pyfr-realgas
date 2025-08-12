<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%inherit file='base'/>
% for mod, name in src_macros:
    <%include file='${mod}'/>
% endfor

<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${fmt(nvars)}]'
              ploc='in fpdtype_t[${fmt(ndims)}]'
              u='in fpdtype_t[${fmt(nvars)}]'
              rcpdjac='in fpdtype_t'>
fpdtype_t src[${nvars}] = {};

% for mod, name in src_macros:
    ${pyfr.expand(name, 't', 'u', 'ploc', 'src')};
% endfor

% for i in range(nvars):
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + src[${i}];
% endfor
</%pyfr:kernel>
