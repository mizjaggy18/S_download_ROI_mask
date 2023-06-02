FROM cytomine/software-python3-base

RUN pip install opencv-python-headless
RUN pip install numpy
RUN pip install shapely


ADD download_roi_mask.py /app/download_roi_mask.py

ENTRYPOINT ["python3", "/app/download_roi_mask.py"]
