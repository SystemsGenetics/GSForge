FROM jupyter/datascience-notebook:latest


#COPY requirements.txt /tmp/
#RUN conda install --yes --file /tmp/requirements.txt && \
#    fix-permissions $CONDA_DIR && \
#    fix-permissions /home/$NB_USER

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install --user $HOME/gsforge&& \
#    fix-permissions $CONDA_DIR && \
#    fix-permissions /home/$NB_USER

COPY ./examples $HOME/examples/

# https://github.com/numba/numba/issues/4032
RUN mkdir /tmp/numba_cache & chmod 777 /tmp/numba_cache & NUMBA_CACHE_DIR=/tmp/numba_cache
