<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%inherit file='base'/>

<%pyfr:kernel name='intconu' ndim='1'
              ulin='in view fpdtype_t[${fmt(nvars)}]'
              urin='in view fpdtype_t[${fmt(nvars)}]'
              ulout='out view fpdtype_t[${fmt(nvars)}]'
              urout='out view fpdtype_t[${fmt(nvars)}]'>
% for i in range(nvars):
% if c['ldg-beta'] == -0.5:
    urout[${i}] = ulin[${i}];
% elif c['ldg-beta'] == 0.5:
    ulout[${i}] = urin[${i}];
% else:
    ulout[${i}] = urout[${i}] = urin[${i}]*${0.5 + c['ldg-beta']}
                              + ulin[${i}]*${0.5 - c['ldg-beta']};
% endif
% endfor
</%pyfr:kernel>
