FROM python:3.11

# Install all required packages
RUN apt install -y \
    ca-certificates \
    curl \
    python3

ENV RFAM_3D=/usr/src/rfam_3d

WORKDIR $RFAM_3D

RUN \
    cd $INFERNAL && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.5.tar.gz && \
    tar -xvzf infernal-1.1.5.tar.gz && \
    cd infernal-1.1.5 && \
    ./configure && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $INFERNAL && \
    rm -r infernal-1.1.5*

RUN curl -sSL https://install.python-poetry.org | python3 -

COPY poetry.lock poetry.lock
COPY pyproject.toml pyproject.toml

RUN PATH="$PATH:/root/.local/bin" poetry config virtualenvs.create false
RUN PATH="$PATH:/root/.local/bin" poetry install

ENTRYPOINT ["/bin/bash"]
