# circle-parameter-estimation
realization of circle parameter estimation based on support vector regression

# Introduction

SVR.m : the main function includes loading datas, constrcuting optimization objectives, constructing constraint condition and processing stage
    
    1 optimization objectives : minimize f(x) = C*H*X
    
    2 constraint condition : A*X <= b
    
    3 processing stage : calculating parameter iteratively based on Augmented Lagrangian Method

# data

There are five group datas in "data" folder, every group contains many point coordinates, and these points 

are located in the edge of a cicle. The algorithm is used to estimate the circle parameter based on the points.
