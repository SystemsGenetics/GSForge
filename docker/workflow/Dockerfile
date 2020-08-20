FROM python:3.8
ENV GSFORGE_INSTALL_CORE=true

RUN groupadd -r gsf_user && useradd --no-log-init -r -g gsf_user gsf_user
USER gsf_user

RUN mkdir $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install --user $HOME/gsforge
