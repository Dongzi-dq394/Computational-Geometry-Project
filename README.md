# Computational-Geometry-Project
This is the repo for my final project of Computational Geometry course, which is mainly about the implmentation and demo of Chan's Algorithm. Chan's Algorithm is an output-sensitive convex hull algorithm that constructs the convex hull for a set of ```n``` points in worst-case optimal ```O(nlogh)``` time and linear space, where ```h``` denotes the number of vertices of convex hull.

* Contributor: [Dongzi Qu](dq394@nyu.edu)
* Instructor: [Prof. Boris Aronov](boris.aronov@nyu.edu)
---
### The description of files:
* chan.py: includes Graham's Scan, Jarvis's March and Binary search for Chan's Algorithm.
* util.py: includes the auxiliary poltting functions for both debugging and demo.

### How to run the algorithm/show the demo:
There are several different ways to run the ```chan.py```. The common style is: ```python chan.py param1 param2 param3 param4```.
* ```param1```: used for testmode, 1 is ```True```, 0 is ```False```. If you want to know how this algorithm work step by step, you should use 1.
* ```param2```: Please choose 2.
* ```param3```: Specifying the size of data set, which is ```n```.
* ```param4```: Set the limit of the number of vertices of the convex hull.

Basically, I will use: ```python chan.py 1 2 15 4```, which means there are 15 points in total and the number of vertices on the boundary will be equal to or smaller than 4. Also, I want to use testmode to show the detail of the running process by several image.
