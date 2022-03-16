FROM python:3.10

ENV RNA /usr/src/rfam
ENV INFERNAL /usr/src/infernal
RUN mkdir $INFERNAL
WORKDIR $RNA

# Install Infernal
RUN \
    cd $INFERNAL && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure --prefix=$INFERNAL/infernal-1.1.2 && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $INFERNAL && \
    rm infernal-1.1.2.tar.gz

ENV PATH="$INFERNAL/infernal-1.1.2/bin:$RNA:$PATH"

ENTRYPOINT ["/bin/bash"]
