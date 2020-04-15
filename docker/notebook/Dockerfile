FROM jupyter/datascience-notebook:latest


#COPY requirements.txt /tmp/
#RUN conda install --yes --file /tmp/requirements.txt && \
#    fix-permissions $CONDA_DIR && \
#    fix-permissions /home/$NB_USER

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install $HOME/gsforge

COPY ./examples $HOME/examples/
