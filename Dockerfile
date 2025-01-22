FROM continuumio/miniconda3

RUN conda update -n base -c defaults -y conda \
&& git clone --branch refactor https://github.com/ozefreitas/M-PARTY.git \
&& conda install -c conda-forge -y mamba \
&& mamba env update -f M-PARTY/ci/ci_environment.yml --name base \
&& bash M-PARTY/ci/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "bin/m-party.py" ]
