FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Install system dependencies for VTK/pyvista
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY . /app

# Install Python dependencies
RUN pip install --no-cache-dir \
    numpy \
    scipy \
    matplotlib \
    pyvista \
    tqdm \
    joblib \
    p-tqdm

# Set default command
CMD ["python", "main.py"]
