FROM python:3.8
ENV GSFORGE_INSTALL_CORE=true

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install --user $HOME/gsforge