FROM jupyter/datascience-notebook:latest

RUN mkdir /opt/gemprospector
COPY . /opt/gemprospector/
RUN pip install /opt/gemprospector/

COPY examples/ examples/
