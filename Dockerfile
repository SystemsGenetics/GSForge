FROM jupyter/datascience-notebook:latest

COPY [^examples]* /opt/gemprospector \
    && pip install /opt/gemprospector

COPY examples/ examples/
