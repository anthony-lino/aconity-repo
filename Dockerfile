# Base image with Python and build tools
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git cmake g++ build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install linear algebra library
RUN apt-get update && apt-get install -y libeigen3-dev cmake

# Set working directory
WORKDIR /opt

# Clone libSLM with submodules
RUN git clone --recursive https://github.com/anthony-lino/libSLM.git

# Install Python bindings
RUN pip install ./libSLM

# Install pyslm
RUN git clone --recursive https://github.com/anthony-lino/pyslm.git 
RUN pip install ./pyslm

# Set the entrypoint to Python by default
ENTRYPOINT ["python3"]
