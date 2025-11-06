FROM python:3.12-slim

ENV IGDATA="/opt/igblast/internal_data"
ENV PATH="/opt/igblast/bin:$PATH"

WORKDIR /app

# Install IgBLAST and system dependencies (including libxml2!)
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    wget \
    curl \
    ca-certificates \
    tar \
    unzip \
    make \
    git \
    libncurses-dev \
    zlib1g \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Download and install IgBLAST
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz && \
    tar -xzf ncbi-igblast-1.22.0-x64-linux.tar.gz && \
    mv ncbi-igblast-1.22.0 /opt/igblast && \
    rm ncbi-igblast-1.22.0-x64-linux.tar.gz

RUN pyir setup

COPY . .

EXPOSE 8000

ENV FLET_NO_RELOAD=1
ENV FLET_SKIP_DAEMON=1
CMD pkill -f "fletd" || true && python Counter/src/main.py