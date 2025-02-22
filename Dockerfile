FROM rocker/r-base:4.3.1
RUN apt-get update && apt-get -y install wget ncbi-blast+
RUN wget https://github.com/ph-u/dNdS_scan/archive/refs/heads/master.zip && unzip master.zip && rm master.zip
RUN mkdir -p data && mkdir -p res
RUN Rscript -e "install.packages('https://CRAN.R-project.org/package=ape&version=5.8', type='source', repos=NULL);install.packages('BiocManager');BiocManager::install('Biostrings', version='3.18')"

##### pipeline settings #####
RUN mv /dNdS_scan-master/src/00_ffn2fa.sh /dNdS_scan-master/binHPC2/
RUN mv /dNdS_scan-master/binHPC2 /
RUN mv /dNdS_scan-master/containerized/* /binHPC2/
RUN for i in `ls /binHPC2/*`;do chmod 755 ${i};done

##### Set env #####
RUN rm -r /dNdS_scan-master
ENV PATH="/binHPC2:${PATH}"
WORKDIR /binHPC2
CMD ["R", "--version"]
# docker rmi -f $(docker images -aq)
# docker image build -t phudock/dnds . && docker run -it phudock/dnds bash
