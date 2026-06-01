# Use explicit lightweight stable base image layer
FROM python:3.11-slim

# Enforce clean working environment workspace path
WORKDIR /app

# Ingest and cache pipeline configurations efficiently
COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Ingest core code elements
COPY app.py app.py

# Expose server entry configuration ports
EXPOSE 5000

# Execute server initialization pass 
CMD ["python", "app.py"]

