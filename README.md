# cEOT
Continuous-time overlap tracker for event data from Neuromorphic Vision Sensors. 

The system takes as input asynchronous events defined by _x_ and _y_ coordinates, and a timestamp, and outputs bounding boxes defined by a center coordinate and its width and height.

# Usage
First, add to path the [MATLAB AER Vision Functions](https://github.com/gorchard/Matlab_AER_vision_functions) from Garrick Orchard.

Then, run the cEOT function:

```
track_cEOT('test_file');
```

The output of the process is an annotation file in text format.