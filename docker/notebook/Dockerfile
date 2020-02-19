# Until an an official 3.8 release comes from below:
#FROM jupyter/datascience-notebook:latest
# A temporary 3.8 image is used from:
# https://gist.github.com/matthewfeickert/9007c41145356897a7173b6d0b9e74f0
# docker run --rm -it -p 10000:8888 systemsgenetics/gsforge_interactive:dev

FROM python:3.8

USER root
WORKDIR /root

SHELL [ "/bin/bash", "-c" ]

RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
        wget \
        curl \
        git \
        make \
        sudo \
        bash-completion && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/*

RUN pip install --upgrade --no-cache-dir pip setuptools wheel && \
    pip install --no-cache-dir numpy jupyter jupyterlab

# Enable tab completion by uncommenting it from /etc/bash.bashrc
# The relevant lines are those below the phrase "enable bash completion in interactive shells"
RUN export SED_RANGE="$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+1)),$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+7))" && \
    sed -i -e "${SED_RANGE}"' s/^#//' /etc/bash.bashrc && \
    unset SED_RANGE

# Create user "docker" with sudo powers
RUN useradd -m docker && \
    usermod -aG sudo docker && \
    echo '%sudo ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers && \
    cp /root/.bashrc /home/docker/ && \
    mkdir /home/docker/data && \
    chown -R --from=root docker /home/docker

# Use C.UTF-8 locale to avoid issues with ASCII encoding
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENV HOME /home/docker

# Have Jupyter notebooks launch without additional command line options
RUN jupyter notebook --generate-config && \
    sed -i -e "/allow_root/ a c.NotebookApp.allow_root = True" ${HOME}/.jupyter/jupyter_notebook_config.py && \
    sed -i -e "/custom_display_url/ a c.NotebookApp.custom_display_url = \'http://localhost:8888\'" ${HOME}/.jupyter/jupyter_notebook_config.py && \
    sed -i -e "/c.NotebookApp.ip/ a c.NotebookApp.ip = '0.0.0.0'" ${HOME}/.jupyter/jupyter_notebook_config.py && \
    sed -i -e "/open_browser/ a c.NotebookApp.open_browser = False" ${HOME}/.jupyter/jupyter_notebook_config.py

RUN chown -R --from=root docker ${HOME}

WORKDIR /home/docker/data
ENV USER docker
USER docker
ENV PATH /home/docker/.local/bin:$PATH


RUN mkdir --parents ${HOME}/gsforge/bin/
COPY . ${HOME}/bin/gsforge/
RUN python -m pip install ${HOME}/bin/gsforge/ --user

# Add the tutorials from GSForge.
COPY ./examples ${HOME}/examples/

# Set the container CMD executable.
CMD [ "/bin/bash" ]