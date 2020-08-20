FROM python:3.7
ENV GSFORGE_INSTALL_CORE=true

#RUN groupadd -r gsf_user && \
# useradd --no-log-init -r -g gsf_user -d /home/gsf_user gsf_user
#USER gsf_user


RUN mkdir -p $HOME/gsforge
COPY . $HOME/gsforge
RUN pip install $HOME/gsforge/setup.py

# https://github.com/numba/numba/issues/4032
RUN mkdir /tmp/numba_cache & chmod 777 /tmp/numba_cache & NUMBA_CACHE_DIR=/tmp/numba_cache
