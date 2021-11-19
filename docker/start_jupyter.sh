#!/bin/bash

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
