FROM python:3.11-slim

WORKDIR /app

COPY pyproject.toml README.md ./
COPY wombat/ wombat/
COPY src/ src/
COPY configs/ configs/

RUN pip install --no-cache-dir .

ENTRYPOINT ["wombat"]
