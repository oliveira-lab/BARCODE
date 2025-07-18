FROM mambaorg/micromamba:1.5.6

# Set the working directory
WORKDIR /workspace

# Copy the environment YAML files and script into the container
COPY environment-pipeline.yml .
COPY environment-analysis.yml .
COPY BARCODE.sh .
COPY models/ /opt/defense-finder-models/

# Create the conda environments
RUN micromamba create -y -n BARCODE_pipeline -f environment-pipeline.yml && micromamba clean --all --yes
RUN micromamba create -y -n BARCODE_analysis -f environment-analysis.yml && micromamba clean --all --yes

# Set the default environment to activate on container start
ENV MAMBA_DEFAULT_ENV=BARCODE_pipeline
ENV PATH=/opt/conda/envs/BARCODE_pipeline/bin:$PATH
ENV DEFENSE_FINDER_MODELS=/opt/defense-finder-models

# Default working directory inside the container
WORKDIR /workspace

# Default command when the container starts
CMD ["/bin/bash"]
