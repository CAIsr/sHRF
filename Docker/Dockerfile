FROM ubuntu:14.04

MAINTAINER Isshaa Aarya <i.aarya@uq.edu.au>

RUN apt-get -qq update && apt-get -qq install -y \
    libpng12-dev libfreetype6-dev libblas-dev liblapack-dev gfortran build-essential xorg \
	curl git unzip make \
	apt-utils 

# make dir for spm installation in /opt/iaSPM

RUN mkdir -p /opt/SPM

WORKDIR /opt/SPM

RUN curl -sSL http://www.fil.ion.ucl.ac.uk/spm/download/restricted/idyll/spm8.zip -o spm8.zip && \
    unzip spm8.zip && \
	rm -rf spm8.zip

ENV PATH "$PATH:/opt/SPM"

# sHRF installation

RUN mkdir -p /opt/sHRF

WORKDIR /opt/sHRF

RUN curl -sSL https://cai.centre.uq.edu.au/files/165/sHRF.zip -o sHRF.zip && \
    unzip sHRF.zip && \
	rm -rf sHRF.zip
	
ENV PATH "$PATH:/opt/sHRF"

# test example .m files
RUN mkdir -p /iaExamples

# set user-defined matlab environment variable
ENV PATH="/usr/local/MATLAB/host/bin:${PATH}"

# set rigid entrypoint for matlab ide
ENTRYPOINT ["matlab"]

# this command is silent at the moment and doesnt do anything
CMD ["-h"]