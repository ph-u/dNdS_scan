FROM rocker/r-base:4.3.1
LABEL org.opencontainers.image.source="https://github.com/ph-u/dNdS_scan"
RUN apt-get update && apt-get -y install wget ncbi-blast+ libxml2-dev
#RUN install.r -e -s -d TRUE BiocManager && install2.r -e -s -d TRUE ape && installBioc.r -e -s -d TRUE Biostrings && rm -rf /tmp/downloaded_packages && rm -rf /var/lib/apt/lists/*
RUN R -q -e "options(warn=2);install.packages('BiocManager', dependencies = T);install.packages('ape', dependencies = T);BiocManager::install('Biostrings', version = 3.18)"
RUN mkdir -p data
#RUN wget 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'

ARG CACHEBUST
RUN echo "$CACHEBUST"
#COPY master.zip /
#RUN unzip master.zip && rm master.zip
RUN wget https://github.com/ph-u/dNdS_scan/archive/refs/heads/master.zip && unzip master.zip && rm master.zip

##### pipeline settings #####
RUN mv /dNdS_scan-master/src/00_ffn2fa.sh /dNdS_scan-master/binHPC2/
#RUN mv /src/00_ffn2fa.sh /binHPC2/
RUN mv /dNdS_scan-master/binHPC2 /
RUN mv /dNdS_scan-master/containerized/* /binHPC2/
#RUN mv /containerized/* /binHPC2/
#RUN mv datasets /binHPC2/
RUN for i in `ls /binHPC2/*`;do chmod 755 ${i};done

##### lib dependency: libgfortran3 #####
# https://gist.github.com/sakethramanujam/faf5b677b6505437dbdd82170ac55322
RUN bash /binHPC2/install-libgfortran3.sh

##### Set env #####
#RUN rm -r /dNdS_scan-master
ENV PATH="/binHPC2:${PATH}"
WORKDIR /binHPC2
CMD ["cp", "/binHPC2/masterTemplate.sh", "masterTemplate.sh"]
