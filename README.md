# rbfinterp2
Local interpolation of data (*x*, *y*, *f*(*x*,*y*)) using polyharmonic spline (PHS) radial basis functions (RBFs).\
Assuming that your data comes from some underlying (possibly unknown) function *f*, the goal of\
*approximation* is to estimate *f*(*x<sub>e*,*y<sub>e*), where (*x<sub>e*, *y<sub>e*) are so-called "evaluation points".  More specifically,\
the goal of *interpolation* (one type of approximation) is to construct a function *g* that matches *f* at every\
node (*x*, *y*), so that *g*(*x<sub>e*,*y<sub>e*) can be used as an estimate for *f*(*x<sub>e*,*y<sub>e*).
## Requirements
* Windows
* Perl
  * Install Strawberry Perl
  * https://strawberryperl.com/
* Python (including numpy, matplotlib, and other packages)
  * Install Anaconda
  * https://www.anaconda.com/
## Instructions
* Download and unzip the repository folder.
  * Push the green "Code" button and select "Download ZIP".
* Open Anaconda Prompt.
  * Press the Start Button in windows, and begin typing "anaconda".  It should pop up.
* Navigate to the folder where you saved this repository.
  * cd C:\Users\gabar\gitRepos\rbfinterp2
### Example problem using random nodes
* Create a new folder for holding coordinates (nodes and evaluation points).
  * mkdir randomCoords
* Create a subfolder of the coordinates folder to store function values at the nodes.
  * mkdir randomCoords\smoothData
* Run three scripts to generate some nodes and function values.
  * perl getNodes.pl randomCoords
  * perl getEvalPts.pl randomCoords
  * perl getFuncVals.pl randomCoords\smoothData
* Run the main script rbfinterp2.pl to interpolate and estimate function values at evaluation points.
  * perl rbfinterp2.pl randomCoords\smoothData y
    * Second input "y" means "yes", the true function IS available for comparison.
### Using your own data
* Create a new folder for holding coordinates (nodes and evaluation points).
  * mkdir coords1
  * Save (delimited by '\n') coordinate data in your coordinate directory.
    * x.txt, y.txt, xe.txt, ye.txt
* Create a subfolder of the coordinates folder to store function values at the nodes.
  * mkdir coords1\data1
  * Save (delimited by '\n') function data in your data directory.
    * f.txt (required), fe.txt (optional, but required for error calculation)
* Run the main script rbfinterp2.pl to interpolate and estimate function values at evaluation points.
  * perl rbfinterp2.pl coords1\data1 y
    * Second input "y" means "yes", the true function IS available for comparison.
