FROM jupyter/datascience-notebook:latest

RUN mkdir /opt/gsforge
COPY . /opt/gsforge/
RUN pip install /opt/gsforge/

COPY ./examples examples/
