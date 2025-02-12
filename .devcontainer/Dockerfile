ARG DEVCONTAINER_BASE=mcr.microsoft.com/devcontainers/base:1.0.9-ubuntu-22.04

#########################################################
# file-normalizer stage
# In order to use BuildKit remote caching, input files must have
# not only the right content hash, but also the right permissions.
# Git only tracks whether the owner can execute a file.
# Here we bring in all files that are going to be used in the
# subsequent stage and normalize the permissions.
#########################################################

FROM --platform=linux/amd64 ${DEVCONTAINER_BASE} as file-normalizer

COPY environment.yml \
     .devcontainer/devcontainer.bashrc  \
     /data/

RUN chmod -R 555 /data/

#########################################################
# devcontainer stage
# Installs all dependencies and tooling for development.
#########################################################

FROM --platform=linux/amd64 ${DEVCONTAINER_BASE} AS devcontainer

# Install needed packages and setup non-root user.
ARG USERNAME="vscode"
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG CONDA_GID=900
ARG CONDA_ENVIRONMENT_NAME=petsird

RUN apt-get update && apt-get install -y \
    libc6-dbg \
    && rm -rf /var/lib/apt/lists/*

# Enable non-root Docker access in container
ARG DOCKER_VERSION="27.0.3"
RUN script=$(curl -fsSL "https://raw.githubusercontent.com/devcontainers/features/2951f0481a488ea43a6f2ea6f18a47f0a0bf744d/src/docker-outside-of-docker/install.sh") \
    && VERSION=${DOCKER_VERSION} DOCKERDASHCOMPOSEVERSION="none" bash -c "$script"

# Setting the ENTRYPOINT to docker-init.sh will configure non-root access to
# the Docker socket if "overrideCommand": false is set in devcontainer.json.
# The script will also execute CMD if you need to alter startup behaviors.
ENTRYPOINT [ "/usr/local/share/docker-init.sh" ]
CMD [ "sleep", "infinity" ]

# Based on https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh -O /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/conda clean --tarballs --index-cache --packages --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && /opt/conda/bin/conda clean --force-pkgs-dirs --all --yes  \
    && groupadd -r conda --gid ${CONDA_GID} \
    && usermod -aG conda ${USERNAME} \
    && chown -R :conda /opt/conda \
    && chmod -R g+w /opt/conda \
    && find /opt -type d | xargs -n 1 chmod g+s

# Create a conda environment from the environment file in the repo root.
COPY --from=file-normalizer --chown=$USER_UID:conda /data/environment.yml /tmp/build/
RUN umask 0002 \
    && /opt/conda/bin/conda env create -f /tmp/build/environment.yml \
    && /opt/conda/bin/conda clean -fy \
    && sudo chown -R :conda /opt/conda/envs

# Add a file that is to be sourced from .bashrc and from the devops pipeline stages
COPY --from=file-normalizer /data/devcontainer.bashrc /opt/devcontainer/

# Add a section to /etc/bash.bashrc that ensures that a section is present at the end of ~/.bashrc.
# We can't just write to .bashrc from here because it will be overwritten if the devcontainer user has
# opted to use their own dotfiles repo. The dotfiles repo is cloned after the postCreateCommand
# in the devcontainer.json file is executed.
RUN echo "\n\
if ! grep -q \"^source /opt/devcontainer/devcontainer.bashrc\" \${HOME}/.bashrc; then\n\
	echo \"source /opt/devcontainer/devcontainer.bashrc\" >> \${HOME}/.bashrc\n\
fi\n" >> /etc/bash.bashrc

ENV CMAKE_GENERATOR=Ninja

# Create a kits file for the VSCode CMake Tools extension, so you are not prompted for which kit to select whenever you open VSCode
RUN mkdir -p /home/vscode/.local/share/CMakeTools \
    && echo '[{"name":"GCC-10","compilers":{"C":"/opt/conda/envs/petsird/bin/x86_64-conda_cos6-linux-gnu-gcc","CXX":"/opt/conda/envs/petsird/bin/x86_64-conda_cos6-linux-gnu-g++"}}]' > /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json \
    && chown vscode:conda /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json

# Install the yardl tool
ARG YARDL_VERSION=0.6.3
RUN wget --quiet "https://github.com/microsoft/yardl/releases/download/v${YARDL_VERSION}/yardl_${YARDL_VERSION}_linux_x86_64.tar.gz" \
    && tar -xzf "yardl_${YARDL_VERSION}_linux_x86_64.tar.gz" \
    && mv yardl "/opt/conda/envs/${CONDA_ENVIRONMENT_NAME}/bin/" \
    && rm "yardl_${YARDL_VERSION}_linux_x86_64.tar.gz"
