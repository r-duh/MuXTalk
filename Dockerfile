FROM python:3.8-slim-buster
WORKDIR /MuXTalk_app
COPY requirements.txt requirements.txt
RUN pip install --upgrade pip
RUN pip3 install -r requirements.txt
COPY . .
ENTRYPOINT ["python3", "run_MuXTalk_v2.py"]