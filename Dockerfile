FROM python:2
MAINTAINER Leighton Pritchard <leighton.pritchard@hutton.ac.uk>

# Install Linux packages and update pip
RUN apt-get update && apt-get install -y \
                               ncbi-blast+ \
			       prodigal \
			       emboss && \
    pip install --upgrade pip && \
    pip install biopython

# Obtain and install primer3 v1.1.4
# For some reason ADD does not give a working .tar.gz file
WORKDIR /primer3
RUN wget https://sourceforge.net/projects/primer3/files/primer3/1.1.4/primer3-1.1.4.tar.gz && \
    tar -xvf primer3-1.1.4.tar.gz
WORKDIR /primer3/src
RUN make all
ENV PATH /primer3/src:${PATH}

# Install find_differential_primers v0.1.3
ADD https://github.com/widdowquinn/find_differential_primers/archive/v0.1.3.tar.gz /fdp/
WORKDIR /fdp
RUN tar -zxvf v0.1.3.tar.gz
WORKDIR /fdp/find_differential_primers-0.1.3
RUN python setup.py install

# Start in host directory
WORKDIR /host_dir

ENTRYPOINT ["find_differential_primers.py"]


