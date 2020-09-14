# Thesis_project
Master's thesis project

Files included in this repository are a part of my thesis project for my master's degree at Faculty of Electrical Engineering and Computing at University of Zagreb. 
I didn't upload the full project but only some of the files to showcase my work. Files included in the C_files folder are a part of a Keil uVision project run on SMT32 
NUCLEO-144 development board and the files included in the Matlab_files folder are some of the MATLAB files used for testing and evaluation of the embedded implementation.


Thesis summary: 

Recognizing Time-Frequency Clustering in Ultrasound Emission Spectra

  Increased water potential in plants conducting tissue, xylem, produces gas bubbles, which when breaking produce mechanical waves of ultrasound frequencies, 
  ultrasound acoustic emissions (UAE). Time-frequency analysis of UAE in real time on an embedded system mounted on a plant and water potential evaluation can 
  greatly improve watering systems in high precision farming. This thesis deals with the development and the implementation on an embedded system of an algorithm 
  set for real time analysis of UAE and transferring compressed data. The algorithm set is made up of OPTICS algorithm, gradient clustering, window merging and 
  result compression. The OPTICS algorithm and gradient clustering are used for extracting clusters. Window merging is used for connecting results of neighboring 
  data blocks, while result compression is used for depicting the recognized clusters with a simple geometric shape. Accuracy of said algorithm set is tested on 
  synthetic data with results over 95 %. The possibility of using the solution on real data is described and methods for adaptive parameter determining have been 
  proposed. Analysis of the implementation on the embedded system SMT32 NUCLEO-144 shows average run time of 994 ms and power consumption of 869 µW.
