FROM jupyter/datascience-notebook:latest

RUN git clone git@github.com:SystemsGenetics/GEMprospector.git /opt/gemprospector \
    && pip install /opt/gemprospector

COPY examples/ examples/
