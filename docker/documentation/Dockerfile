FROM jupyter/datascience-notebook:latest
RUN Rscript -e 'install.packages("BiocManager", dependencies=TRUE, repos="http://cran.rstudio.com/")' &&\
    Rscript -e 'BiocManager::install(c("edgeR", "DESeq2", ask = F))'


COPY . $HOME/gsforge
ENV PATH="$HOME$USER/.local/bin:${PATH}"
RUN pip install --user $HOME/gsforge[docs,R_support]
RUN ipython kernel install --name "gsfenv" --user

# https://github.com/numba/numba/issues/4032
RUN mkdir /tmp/numba_cache & chmod 777 /tmp/numba_cache & NUMBA_CACHE_DIR=/tmp/numba_cache

CMD rm -rf ${GSFORGE_DEMO_DATA}/builtdocs/ &&\
    mkdir -p $HOME/gsforge/docs/source/API/ &&\
    sphinx-apidoc --separate --force -o docs/source/API/ $HOME/gsforge &&\
    python -m sphinx ./gsforge/docs/source/ ${HOME}/builtdocs/ &&\
    touch ${HOME}/builtdocs/.nojekyll &&\
    cp -rf ${HOME}/builtdocs/ ${GSFORGE_DEMO_DATA}/builtdocs/
