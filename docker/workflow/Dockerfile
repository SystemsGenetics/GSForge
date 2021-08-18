FROM python:3.7
ENV GSFORGE_INSTALL_CORE=true

#RUN groupadd -r gsf_user && \
# useradd --no-log-init -r -g gsf_user -d /home/gsf_user gsf_user
#USER gsf_user


RUN mkdir /opt/gsforge
COPY . /opt/gsforge/
RUN pip install /opt/gsforge/

# https://github.com/numba/numba/issues/4032
RUN mkdir /tmp/numba_cache & chmod 777 /tmp/numba_cache & NUMBA_CACHE_DIR=/tmp/numba_cache
