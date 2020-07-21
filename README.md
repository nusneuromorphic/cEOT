# cEOT
Continuous-time overlap tracker for event data from Neuromorphic Vision Sensors. 

<p align="center">
  <a href="https://drive.google.com/file/d/1gRb1eC2RDM0ZMFhPZQ2mFYq_AulbJXzj/preview">
    <img src="https://i.ibb.co/qd33x5m/Slide2-PB.jpg" alt="cEOT Video Demo" width="600"/>
  </a>
</p>

The system takes as input asynchronous events defined by _x_ and _y_ coordinates, and a timestamp, and outputs bounding boxes defined by a center coordinate and its width and height.

# Usage
First, add to path the [MATLAB AER Vision Functions](https://github.com/gorchard/Matlab_AER_vision_functions) from Garrick Orchard.

Then, run the cEOT function:

```
track_cEOT('test_file', 0, 0);
```

The output of the process is an annotation file in text format.

An example script is included along with a short sample recording in bin format. 
