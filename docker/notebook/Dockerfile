FROM jupyter/datascience-notebook:latest

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install --user $HOME/gsforge


#COPY ./examples $HOME/examples/

# https://github.com/numba/numba/issues/4032
RUN mkdir /tmp/numba_cache & chmod 777 /tmp/numba_cache & NUMBA_CACHE_DIR=/tmp/numba_cache
