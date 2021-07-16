.. ****************************** Begin license ********************************
.. Copyright (C) Laboratory of Imaging technologies,
..               Faculty of Electrical Engineering,
..               University of Ljubljana.
..
.. This file is part of PyXOpto.
..
.. PyXOpto is free software: you can redistribute it and/or modify
.. it under the terms of the GNU General Public License as published by
.. the Free Software Foundation, either version 3 of the License, or
.. (at your option) any later version.
..
.. PyXOpto is distributed in the hope that it will be useful,
.. but WITHOUT ANY WARRANTY; without even the implied warranty of
.. MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
.. GNU General Public License for more details.
..
.. You should have received a copy of the GNU General Public License
.. along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
.. ******************************* End license *********************************

PyXOpto documentation
=====================

This is `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_-based
documentation. Follow these steps to build the documentation from source.

Full build / rebuild of the HTML documentation
----------------------------------------------

First, render the documentation templates by running the following command
from the root documentation directory (the one that includes the conf.py
and render.py files):

.. code:: bash

    python render.py

The next step will process the python source files of the PyXOpto project.

Linux OS
^^^^^^^^

Run :code:`custom_build.sh` from the root directory of documentation
(the one with the :code:`custom_build.sh` file).

.. code-block:: bash

    ./custom_build.sh

Windows OS
^^^^^^^^^^

Run :code:`custom_build.bat` from the root directory of documentation
(the one with the :code:`custom_build.bat` file).

.. code-block:: bash

    custom_build.bat

This will build / extract the documentation from the Python source code and
create a number of .rst files in the :code:`apidoc` directory and finally build
the HTML documentation.

Rebuilding after changing documents in the :code:`source` directory
-------------------------------------------------------------------

When rebuilding the documentation due to changes of documents in the
:code:`source` directory, it is sufficient to run :code:`make html`.

.. code-block:: bash

    make html

This will take the existing files from the :code:`apidoc` directory and updated
the HTML in a fraction of time required by the full build. Note that any
changes to the documentation in the source code of PyXOpto will not be visible
in the HTML.

If any of the documentation template files listed in :py:mod:`render.py` are
changed, the templates need to be reprocessed and the HTML rebuild:

.. code:: bash

    python render.py
    make html
