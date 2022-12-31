# rbfinterp2
Local interpolation of data (x, y, f(x,y)) using polyharmonic spline (PHS) radial basis functions (RBFs).\
Assuming that your data comes from some underlying (possibly unknown) function f, the goal of\
approximation is to estimate f(xe,ye), where (xe, ye) are so-called "evaluation points".  More specifically,\
the goal of interpolation (one type of approximation) is to construct a function g that matches f at every\
node (x, y), then use g(xe,ye) to estimate f(xe,ye).
## Requirements
* Windows
* Perl
  * Install Strawberry Perl
  * https://strawberryperl.com/
* Python
  * Install Anaconda
  * https://www.anaconda.com/
## Instructions
After each step in the instructions, I will try to give an example.
* Download and unzip the repository folder
* Open Anaconda Prompt
  * Press the Start button in windows, and begin typing "anaconda".  It should pop up.
* Navigate to the folder where you saved this repo
  * cd C:\Users\gabar\gitRepos\rbfinterp2
* Create a new folder for holding coordinates
  * mkdir randomCoords
* Create a subfolder of the coordinates folder for holding function values at the coordinates.
  * mkdir randomCoords\smoothData
* Run three scripts to generate some nodes and function values.
  * perl getRandomNodes.pl randomCoords
  * perl getEvalPts.pl randomCoords
  * perl getSmoothData.pl randomCoords\smoothData
* Run the main script rbfinterp2.pl to interpolate and estimate function values at evaluation points.
  * perl rbfinterp2.pl randomCoords\smoothData y
    * Second input "y" means "yes", the true function IS available for comparison
