# Track Model Optimization

The primary task of the track model is to provide the simulated vehicle with a path to follow. It is generally done by defining the corner radius or track curvature as a function of distance. <br>
The following track optimization tool takes the raw GPS coordinates (lattitude, longitude, altitude) of points along the left and right boundaries (according to the driving direction) and converts them into a ribbon model for the track in frenet coordinates. The ribbon model is passed through an optimizer which allows for smoother track model with continuous curvature, angles and track widths which can be used with the vehicle and tire models.
