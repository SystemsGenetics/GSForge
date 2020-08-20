FROM jupyter/datascience-notebook:latest


#COPY requirements.txt /tmp/
#RUN conda install --yes --file /tmp/requirements.txt && \
#    fix-permissions $CONDA_DIR && \
#    fix-permissions /home/$NB_USER

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install --user $HOME/gsforge&& \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

COPY ./examples $HOME/examples/
