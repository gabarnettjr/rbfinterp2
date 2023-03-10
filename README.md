# rbfinterp2
Local interpolation of data (*x*, *y*, *f*(*x*,*y*)) using polyharmonic spline (PHS) radial basis functions (RBFs).\
Assuming that your data comes from some underlying (possibly unknown) function *f*, the goal of\
*approximation* is to estimate *f*(*x<sub>e*,*y<sub>e*), where (*x<sub>e*, *y<sub>e*) are so-called "evaluation points".  More specifically,\
the goal of *interpolation* (one type of approximation) is to construct a function *g* that matches *f* at every\
node (*x*, *y*), so that *g*(*x<sub>e*,*y<sub>e*) can be used to estimate *f*(*x<sub>e*,*y<sub>e*).
## Requirements
* Python (including numpy, matplotlib, and other packages)
  * https://www.anaconda.com/
## Instructions
* Download and unzip the repository folder.
  * Push the green "Code" button and select "Download ZIP".
* Open Anaconda Prompt.
  * In windows press the Start Button and begin typing "anaconda".  It should pop up.
  * In mac or linux, open it like you normally open stuff.
* NOTE: In the instructions that follow, if you are using mac or linux, then\
  replace each backslash \\ with a forward slash /.
* Navigate to the folder where you saved this repository.
  * cd C:\Users\gabar\gitRepos\rbfinterp2
### Using random nodes and a smooth underlying function (toy problem)
* Create a new folder for holding coordinates (nodes and evaluation points).
  * mkdir randomCoords
* Create a subfolder of the coordinates folder to store function values at the nodes.
  * mkdir randomCoords\smoothData
* Run three scripts to generate some nodes and function values.
  * python getNodes.py randomCoords
  * python getEvalPts.py randomCoords
  * python getFuncVals.py randomCoords\smoothData
* Run the main script rbfinterp2.py to interpolate and estimate function values at evaluation points.
  * python rbfinterp2.py randomCoords\smoothData y
    * Second input "y" means "yes", the true function IS available for comparison.
### Using your own data (real problem)
* Create a new folder for holding coordinates (nodes and evaluation points).
  * mkdir coords1
  * Save (delimited/separated by newlines) coordinate data in your coordinate directory.
    * x.txt, y.txt, xe.txt, ye.txt
* Create a subfolder of the coordinates folder to store function values at the nodes.
  * mkdir coords1\data1
  * Save (delimited/separated by newlines) function data in your data directory.
    * f.txt (required), fe.txt (optional, but required for error calculation)
* Run the main script rbfinterp2.py to interpolate and estimate function values at evaluation points.
  * python rbfinterp2.py coords1\data1 n
    * Second input "n" means "no", the true function is NOT available for comparison.
## More Help
Navigate to where you saved the repo and execute this command.
* python rbfinterp2.py --help
