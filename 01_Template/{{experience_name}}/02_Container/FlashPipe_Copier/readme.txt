








docker build -t flashpipe_copier .


docker run --name flashpipe_copier --rm -it -v /mnt:/mnt flashpipe_copier /bin/bash


