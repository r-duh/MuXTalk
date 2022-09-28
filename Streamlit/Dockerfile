FROM python:3.8-slim-buster
EXPOSE 8501
WORKDIR /MuXTalk_Streamlit_app
RUN apt-get update && apt-get install -y \
    build-essential \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*
COPY . .
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt
ENTRYPOINT ["streamlit", "run", "MuXTalk_viz_Streamlit_Docker.py", "--server.port=8501", "--server.address=0.0.0.0"]