from __future__ import print_function, unicode_literals, absolute_import, division


import sys
import numpy as np
import os
import cytomine
# import matplotlib.pyplot as plt
import re
# import pandas as pd
# import skimage.exposure
from shapely.geometry import Polygon, Point, mapping, shape, box
from shapely.affinity import affine_transform
# from shapely import wkt
# from glob import glob
# from tifffile import imread
# from cytomine import Cytomine, models, CytomineJob
from cytomine.models import AbstractImage, Annotation, AnnotationTerm, AnnotationCollection, ImageInstanceCollection, Job, User, JobData, Project, ImageInstance, Property
from cytomine.models.ontology import Ontology, OntologyCollection, Term, RelationTerm, TermCollection

# from csbdeep.utils import Path, normalize
# from PIL import Image

# import matplotlib.pyplot as plt
# import time
import cv2
# import math
# import csv

from argparse import ArgumentParser
# import json
import logging
import logging.handlers
import shutil

__author__ = "WSH Munirah W Ahmad <wshmunirah@gmail.com>"
__version__ = "1.0.0"
# Download ROI mask from nuclei annotations (Date created: 01 June 2023)

def prep_coor(coord):
    regex="[-+]?(?:\d*\.\d+|\d+)"
    coord2=re.findall(regex,coord)
    # coord2=int(str(coord2))
    # print('coord2: ',coord2)
    # print(type(coord2))

    x=[]
    y=[]
    for i, coor in enumerate(coord2):
        coor=float(coor)
        # print(type(coor))
        if (i % 2)==0:
            x.append(coor)
        if (i % 2)==1:
            y.append(coor)
    # print(x)
    # print(y)
    return x,y, coord2
    # return coord2

def process_nuclei(mask,r1,c1,polygon,label,mask_only):
        # print(roi_height)
                    
    point2=polygon
    x, y = point2.exterior.coords.xy
    # print(point2)
    # print(x)
    x=[int(k)-int(c1) for k in x]
    y=[int(r1)-int(k) for k in y]

    # print(x)
    yy=[]
    for k in range(0,len(x)-1):
        yy.append([float(x[k]),float(y[k])])
    # print(yy)

    point3=Polygon(yy)
    # point3=affine_transform(point3, [1, 0, 0, -1, 0, roi_height])
    # f.write("{};{};{};{};{};{};{};{}\n".format(annotation2.id,annotation2.image,annotation2.project,annotation2.term,annotation2.user,annotation2.area,annotation2.perimeter,str(point3)))

    int_coords = lambda x: np.array(x).round().astype(np.int32)
    # print(int_coords)
    exterior = [int_coords(point3.exterior.coords)]
    # print(exterior)   

    isClosed = True

    # Blue color in BGR
    if mask_only==1:
        color = (0,0,0)
    else:
        color = (255, 0, 0)

    # Line thickness of 3 px
    thickness = 3
    # Using cv2.polylines() method
    # Draw a Blue polygon with 
    # thickness of 1 px
    image_ = cv2.polylines(mask, exterior, 
                        isClosed, color, thickness)                    

    # filled = np.zeros_like(roi_mask)
    filled = cv2.fillPoly(image_, exterior, color=label)
    # print(np.amax(filled))
    # result = skimage.exposure.rescale_intensity(filled, in_range=(127.5,255), out_range=(0,255)).astype(np.uint8)
    return filled     

def run(cyto_job, parameters):
    logging.info("----- Download ROI Mask v%s -----", __version__)
    logging.info("Entering run(cyto_job=%s, parameters=%s)", cyto_job, parameters)

    job = cyto_job.job
    user = job.userJob
    project = cyto_job.project
    mask_only = parameters.mask_only
    reviewed_only = parameters.reviewed_only
    id_roi_term=parameters.id_roi_term
    id_cell_term=parameters.id_cell_term
    id_c0_term=parameters.id_c0_term
    id_c1_term=parameters.id_c1_term
    id_c2_term=parameters.id_c2_term
    id_c3_term=parameters.id_c3_term

    terms = TermCollection().fetch_with_filter("project", parameters.cytomine_id_project)
    job.update(status=Job.RUNNING, progress=1, statusComment="Terms collected...")
    print(terms)

    #Select images to process
    images = ImageInstanceCollection().fetch_with_filter("project", project.id)       
    list_imgs = []
    if parameters.cytomine_id_images == 'all':
        for image in images:
            list_imgs.append(int(image.id))
    else:
        list_imgs = parameters.cytomine_id_images
        list_imgs2 = list_imgs.split(',')
        
    print('Print list images:', list_imgs2)

    working_path = os.path.join("tmp", str(job.id))
   
    if not os.path.exists(working_path):
        logging.info("Creating working directory: %s", working_path)
        os.makedirs(working_path)
    try:

        # id_project=project.id   
        # output_path = os.path.join(working_path, "mask_images.zip")
                
        # output_path = os.path.join(working_path, "classification_results.csv")
        # f= open(output_path,"w+")

        # f.write("AnnotationID;ImageID;ProjectID;JobID;TermID;UserID;Area;Perimeter;Hue;Value;WKT \n")
        
        #Go over images
        for id_image in list_imgs2:

            print('Current image:', id_image)
            imageinfo=ImageInstance(id=id_image,project=parameters.cytomine_id_project)
            imageinfo.fetch()
            # im_width=imageinfo.width
            # im_height=imageinfo.height

            annotation_params = {
                "project": parameters.cytomine_id_project,
                "image": id_image,
                "showWKT": True,
                "showMeta": True,
                "showTerm": True,
                "showImage": True,
                "showGIS": True,                
                "annotation": True,
                
            }
            
            if reviewed_only==1:
                review_annotations = AnnotationCollection(**annotation_params, reviewed=True).fetch()
                annotations = review_annotations   
            else:
                user_annotations = AnnotationCollection(**annotation_params).fetch()
                algo_annotations = AnnotationCollection(**annotation_params, includeAlgo=True).fetch()
                annotations = user_annotations + algo_annotations

            logging.debug(annotations)
            print(annotations)

            roi_num=0
            roi_rc={}
            roi_poly_all={}

            # roi_path=os.path.join(working_path,str(annotations.project)+'/')
            # f= open(roi_path+"/"+str(id_image)+".csv","w+")
            # f.write("ID;Image;Project;Term;User;Area;Perimeter;WKT \n")

            for annotation in annotations:
                # image=annotation.image
                # print(type(image))
                term=annotation.term
                # term=str(term)
                # regex = '\d+'          
                # term = re.findall(regex, term) 
                # print(term)

                # print(type(term),type([id_roi_term])) 
                # print(term,[id_roi_term])
                
                # refine=[]
                if term==[id_roi_term]: #str([1486]): #ROI term [1486]; ROI-WSI term 145085; Nuclei_c1weak 114720; nuclei 1517
                    roi_num=roi_num+1
                    print('roi num: ',roi_num)

                    # f.write("{};{};{};{};{};{};{};{}\n".format(annotation.id,annotation.image,annotation.project,annotation.term,annotation.user,annotation.area,annotation.perimeter,annotation.location))

                    # roi_png_filename=os.path.join(roi_path+str(annotation.id)+'.png')
                    # annotation.dump(dest_pattern=roi_png_filename)
                    # im = cv2.imread(roi_png_filename)

                    coord=annotation.location                
                    # print('coord: ',coord)
                    # print(type(coord))

                    x,y,coord2=prep_coor(coord)
                    roi_width=int(x[1]-x[0])
                    # print(roi_width)
                    roi_height=int(y[1]-y[2])
                    # print(roi_height)
                    r1=y[0]
                    c1=x[0]
                    roi_rc[roi_num]=r1,c1
                    print('type roi_rc: ',type(roi_rc))
                    print('roi_rc: ',roi_rc)
                    
                    roi_mask=np.zeros((roi_height,roi_width))
                    mask=np.zeros((roi_height,roi_width,3))
                    
                    roi_poly=Polygon([(x[0],y[0]),(x[1],y[1]),(x[2],y[2]),(x[3],y[3]),(x[4],y[4])])
                    print('ROI Poly: ',roi_poly)
                    roi_poly_all[roi_num]=roi_poly
                    # print(roi_poly_all)

                    filled=[]

                    for i, annotation2 in enumerate(annotations):
                        term=annotation2.term
                        # term=str(term)
                        # regex = '\d+'          
                        # term = re.findall(regex, term) 
                        # print(term)

                        coord=annotation2.location
                        # print('coord nuclei: ',coord)
                        x,y,coord2=prep_coor(coord)
                        yy=[]
                        for k in range(0,len(coord2)-1,2):
                            yy.append([float(coord2[k]),float(coord2[k+1])])

                        polygon=Polygon(yy)
                        # print('polygon: ',polygon)

                        if roi_poly.contains(polygon):
                                    
                            # nucleipoints=[]                            
                            
                            if mask_only==1:
                                if term==[id_cell_term] or term==[id_c0_term] or term==[id_c1_term] or term==[id_c2_term] or term==[id_c3_term]:
                                    label = (255, 255, 255) 
                                    filled=process_nuclei(mask,r1,c1,polygon,label,mask_only)    

                            else:

                                #BGR label (N:blue,W:green,M:yellow,S:red)                            
                                if term==[id_c0_term]:#Nuclei_c0negative 114743
                                    # print(i)
                                    label = (255, 159, 0) 
                                    filled=process_nuclei(mask,r1,c1,polygon,label,mask_only)

                                if term==[id_c1_term]:#Nuclei_c1weak 114720
                                    # print(i)
                                    label = (0, 255, 0) 
                                    filled=process_nuclei(mask,r1,c1,polygon,label,mask_only)
                                
                                if term==[id_c2_term]:#Nuclei_c2moderate 114751
                                    # print(i)
                                    label = (0, 216, 255) 
                                    filled=process_nuclei(mask,r1,c1,polygon,label,mask_only)
                                
                                if term==[id_c3_term]:#Nuclei_c3strong 114759
                                    # print(i)
                                    label = (0, 0, 255) 
                                    filled=process_nuclei(mask,r1,c1,polygon,label,mask_only)
                               
                    mask_filename=os.path.join(working_path,str(annotation.id)+'_mask.png')
                    if len(filled)!=0:
                        # print(filled)
                        # cv2.imshow('image',filled)
                        # cv2.waitKey(0)
                        print(mask_filename)
                        cv2.imwrite(mask_filename,filled)



        shutil.make_archive("mask_images", 'zip', working_path)
        job_data = JobData(job.id, "Generated File", "mask_images.zip").save()
        job_data.upload("mask_images.zip")
        job.update(progress=100, statusComment="Masks saved!")


    finally:
        logging.info("Deleting folder %s", working_path)
        shutil.rmtree(working_path, ignore_errors=True)
        logging.debug("Leaving run()")


    job.update(status=Job.TERMINATED, progress=100, statusComment="Finished.") 

if __name__ == "__main__":
    logging.debug("Command: %s", sys.argv)

    with cytomine.CytomineJob.from_cli(sys.argv) as cyto_job:
        run(cyto_job, cyto_job.parameters)