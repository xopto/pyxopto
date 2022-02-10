#!/bin/bash
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

# (re)build the RNG support functions to slightly reduce the MC overhead
python3 -c "import xopto; xopto.rebuild(verbose=False)"

if [ $EUID -eq 0 ]; then
	if [ -e "/dev/dri/renderD128" ]; then
		egrep -i "^render" /etc/group;
		if [ $? -eq 0 ]; then
			echo "Render group exists"
		else
			echo "Render group does not exist ... creating"
			if [ "$(stat -c "%-G" /dev/dri/renderD128)" == "UNKNOWN" ]
			then
				gid=$(stat -c "%-g" /dev/dri/renderD128)

				echo "Creating group \"render\" with id=$gid"
				groupadd -g $gid render

				echo "Adding render group to user \"$NB_USER\""
				usermod -a -G render $NB_USER
				echo "Groups of user \"$NB_USER\":"
				groups $NB_USER
			fi
		fi
	fi
	su $NB_USER -c "jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser"

else
	jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
fi
