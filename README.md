# rbfinterp2
Local interpolation of data (x, y, f(x,y)) using polyharmonic spline (PHS) radial basis functions (RBFs).\
Assuming that your data comes from some underlying (possibly unknown) function f, the goal of\
approximation is to estimate f(xe,ye), where (xe, ye) are so-called "evaluation points".  More specifically,\
the goal of interpolation (one type of approximation) is to construct a function g that matches f at every\
node, then estimate f(xe,ye) by evaluating g(xe,ye).
## Requirements
* Windows
* Perl
  * Install Strawberry Perl
  * https://strawberryperl.com/
* Python
  * Install Anaconda
  * https://www.anaconda.com/
## Instructions
Here are some instructions for how to test out this repo.  After some of the instructions, I will give an example.
* Download and unzip the repository folder
* Open Anaconda Prompt
  * Press the start button and start typing "anaconda".  It should pop up.
* Navigate to the folder where you saved this repo
  * cd C:\Users\gabar\gitRepos\phsinterp2
* Create a new folder for holding coordinates
  * mkdir randomCoords
* Create a subfolder of the coordinates folder for holding function values at the coordinates.
  * mkdir randomCoords\smoothData
* Run the script getRandomNodes.pl to generate some "random" nodes.
  * perl getRandomNodes.pl randomCoords
* Run the script getEvalPts.pl to generate some evaluation points.
  * perl getEvalPts.pl randomCoords
* Run the script getSmoothData.pl to sample a smooth function on your nodes.
  * perl getSmoothData.pl randomCoords\smoothData
* Run the script rbfinterp2.pl to interpolate and estimate function values at the evaluation points.
  * perl rbfinterp2.pl randomCoords\smoothData y
