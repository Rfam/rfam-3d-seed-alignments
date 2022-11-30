FROM python:3.10

ENV RNA /usr/src/rfam
ENV INFERNAL /usr/src/infernal
ENV SCRIPTS /usr/src/scripts
RUN mkdir $INFERNAL
RUN mkdir $SCRIPTS
WORKDIR $RNA

# Install Infernal
RUN \
    cd $INFERNAL && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.4.tar.gz && \
    tar -xvzf infernal-1.1.4.tar.gz && \
    cd infernal-1.1.4 && \
    ./configure --prefix=$INFERNAL/infernal-1.1.4 && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $INFERNAL && \
    rm infernal-1.1.4.tar.gz

# Install reformatting scripts
RUN cd $SCRIPTS && git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    chmod +x jiffy-infernal-hmmer-scripts/ali-pfam-lowercase-rf-gap-columns.pl

ADD requirements.txt .

RUN pip install -r requirements.txt

ENV PATH="$INFERNAL/infernal-1.1.4/bin:$RNA:$PATH"
ENV PATH="$SCRIPTS/jiffy-infernal-hmmer-scripts:$RNA:$PATH"

ENTRYPOINT ["/bin/bash"]
