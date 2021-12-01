# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################
import sys

from typing import Tuple

import jinja2

def render(templates: Tuple[str, str], context: dict, verbose: bool = False):
    '''
    Renders template files with given context:

    Parameters
    ----------
    templates: Tuple[str, str]
        A pair of strings that define the source template file and the
        destination file for the rendered template.
    context: dict
        Context that is used to render the template file.
    verbose: bool
        Turns on verbose reporting to stdout.
    '''
    for src, dest in templates:
        src_content = None
        with open(src) as fid:
            if verbose:
                print('Loading template from:', src)
            src_content = fid.read()
        if src_content is not None:
            template = jinja2.Template(src_content)
            out_content = template.render(context)
            with open(dest, 'w') as fid:
                if verbose:
                    print('Saving rendered template to:', dest)
                fid.write(out_content)

MCML_CONTEXT = {'MC':'mcml'}
''' Context used to render the mcml documentation tempates '''
MCML_TEMPLATES = [
    ('source/mcbase/trace.template.rst', 'source/mcml/trace.rst'),
    ('source/mcbase/fluence.template.rst', 'source/mcml/fluence.rst'),
    ('source/mcbase/phase_function.template.rst', 'source/mcml/phase_function.rst'),
    ('source/mcbase/sampling_volume.template.rst', 'source/mcml/sampling_volume.rst'),
    ('source/mcbase/simulator_data_types.template.rst', 'source/mcml/simulator_data_types.rst'),
    ('source/mcbase/simulator_options.template.rst', 'source/mcml/simulator_options.rst'),
]
''' Mcml documentation templates and corresponding output files '''

MCVOX_CONTEXT = {'MC':'mcvox'}
''' Context used to render the mcvox documentation tempates '''
MCVOX_TEMPLATES = [
    ('source/mcbase/trace.template.rst', 'source/mcvox/trace.rst'),
    ('source/mcbase/fluence.template.rst', 'source/mcvox/fluence.rst'),
    ('source/mcbase/phase_function.template.rst', 'source/mcvox/phase_function.rst'),
    ('source/mcbase/sampling_volume.template.rst', 'source/mcvox/sampling_volume.rst'),
    ('source/mcbase/simulator_data_types.template.rst', 'source/mcvox/simulator_data_types.rst'),
    ('source/mcbase/simulator_options.template.rst', 'source/mcvox/simulator_options.rst'),

]
''' Mcvox documentation templates and corresponding output files '''

MCCYL_CONTEXT = {'MC':'mccyl'}
''' Context used to render the mccyl documentation tempates '''
MCCYL_TEMPLATES = [
    ('source/mcbase/trace.template.rst', 'source/mccyl/trace.rst'),
    ('source/mcbase/fluence.template.rst', 'source/mccyl/fluence.rst'),
    ('source/mcbase/phase_function.template.rst', 'source/mccyl/phase_function.rst'),
    ('source/mcbase/sampling_volume.template.rst', 'source/mccyl/sampling_volume.rst'),
    ('source/mcbase/simulator_data_types.template.rst', 'source/mccyl/simulator_data_types.rst'),
    ('source/mcbase/simulator_options.template.rst', 'source/mccyl/simulator_options.rst'),

]
''' Mccyl documentation templates and corresponding output files '''

if __name__ == '__main__':
    render(MCML_TEMPLATES, MCML_CONTEXT, verbose=True)
    render(MCCYL_TEMPLATES, MCCYL_CONTEXT, verbose=True)
    render(MCVOX_TEMPLATES, MCVOX_CONTEXT, verbose=True)
