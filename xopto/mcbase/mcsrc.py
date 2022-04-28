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

from typing import List
import re
import os.path
import jinja2
import time


RE_LICENSE = re.compile(
    '/\*+ Begin license '
    '[ a-zA-Z0-3.,;:,.\*\n\r\t()/<>]* End license \*+/')
''' Compiled regular expression that matches the file header - license text. '''


def fuse(inputs: List[str], output: str, striplicense: bool = False,
         verbose: bool = False) -> str:
    '''
    Fuses the source files in inputs into the template defined by the output.

    Parameters
    ----------
    inputs: List[str]
        A list of source files to fuse.
    output: str
        Output template file.
    striplicense: bool
        Strip the license headers if set to True.
    verbose: bool
        Turns on verbose reporting.

    Note
    ----
    The include directives are excluded/commented before fusing the content
    of source files into the template.
    The data context keys are obtained by replacing all dots in the source file
    names with underscores. The output template file should therefore reference
    the source files as "src_h" instead of "src.h".
    '''
    context = {}

    t_start = time.perf_counter()

    def replace_inlude(match):
        # print('replacing', match)
        return '/* {} */'.format(match.group(0))

    for fullpath in inputs:
        src = ''
        filename = os.path.basename(fullpath)
        with open(fullpath, 'r') as fid:
            src = fid.read()

        if striplicense:
            src = RE_LICENSE.sub('', src)

        clean_src, _ = re.subn(
            r'[ \t]*#include[ \t]+"[_a-zA-Z.]+"[ \t]*', replace_inlude, src)


        var = filename.replace('.', '_')
        var = var.replace('-', '_')
        context[var] = clean_src

    with open(output, 'r') as fid:
        template = jinja2.Template(fid.read())

    result = template.render(context)

    if verbose:
        print('OpenCL source files fused in {:3f} ms.'.format(
            (time.perf_counter() - t_start)*1e3))

    return result
    

if __name__ == '__main__':
    from xopto import MCBASE_PATH, MCML_PATH

    _MCML_FUSION_SPEC = {
        'inputs': [
            os.path.join(MCBASE_PATH, 'kernel', 'mcbase.template.h'),
            os.path.join(MCBASE_PATH, 'kernel', 'mcbase.template.c'),
            os.path.join(MCML_PATH, 'kernel', 'mcml.template.h'),
            os.path.join(MCML_PATH, 'kernel', 'mcml.template.c'),
        ],
        'output': os.path.join(MCML_PATH, 'kernel', 'mcml.fusion.template.h')
    }

    r = fuse(
        **_MCML_FUSION_SPEC    
    )
    
    from xopto import USER_TMP_PATH
    with open(os.path.join(USER_TMP_PATH, 'fused.h'), 'w') as fid:
        fid.write(r)
