# Base image with Python and build tools
FROM python:3.9

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gdb \
    vim \
    git cmake g++ build-essential \
    python3-dbg \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN ulimit -c unlimited

ENV COREFILE=/tmp/core.%e.%p

# Install linear algebra library
RUN apt-get update && apt-get install -y libeigen3-dev cmake

# Install Blender CLI (headless)
RUN apt-get update && apt-get install -y \
    blender \
    xvfb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /repos

# Clone libSLM with submodules
RUN git clone --recursive https://github.com/anthony-lino/libSLM.git

# Install Python bindings
RUN pip install ./libSLM

# Install pyslm
RUN git clone --recursive https://github.com/anthony-lino/pyslm.git 
RUN pip install ./pyslm

# Set the entrypoint to Python by default
ENTRYPOINT ["python3"]
