
This image contains:

 - R 4.4
 - Rstudio server

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t qc_r .

# ######################
     RUN THE IMAGE
# ######################

docker run -d --name qc_r -p 8787:8787 -v /mnt:/mnt -e DEFAULT_USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -e PASSWORD=yourpass qc_r


