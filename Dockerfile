FROM cytomine/software-python3-base

RUN pip install numpy
RUN pip install shapely
RUN pip install cv2

ADD download_roi_mask.py /app/download_roi_mask.py

ENTRYPOINT ["python3", "/app/download_roi_mask.py"]