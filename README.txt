
Rocket simulator created by John Nguyen

Requirements:

Python 3 with VPython

Usage:

	Launch vehicles are defined by a FreeRocket object, which, besides being initialized with numerous properties including mass, drag coefficient, and moments of inertia, also require a thrust function, a FinSet object, and a WindProfile object.  Once all these are provided and the FreeRocket successfully initialized, all the user must do is run FreeRocket.simulate().  This performs all numerical integration and graphing functions in a self-contained loop.  Then the user may run FreeRocket.flightReport() to print summary values such as apogee, maximum speed, max-Q, and ground hit velocity.  Several examples are provided near the beginning of the simulation loop.

	One must note that "FreeRocket" refers to the fact that they are aerodynamically stabilized with no guidance, not an association with free and open-source software.  Though of course, as this is an MIT-licensed program, the latter interpretation is not completely without merit.

	The simulation loop can be rewritten however the user wishes to.  One common optimization is to use a very fine time step (<0.01s) until the rocket clears the launch rail, then switch to a coarser time step (~0.05s) for the powered ascent, then ~0.1s for the unpowered ascent, then 0.5s or more for the recovery portion.  Simply call the simulate() method of each FreeRocket that you wish to step within a loop.

	Multiple FreeRocket objects may be instantiated and simulated at once.  Currently, due to the lack of guidance and control, this is only really useful for multi-stage sounding rockets or perhaps a "drag-race" style simultaneous launch.  However, when guidance is implemented, a wide variety of use-cases will open up for multi-vehicle simulation, such as rendezvous and docking, interception, and more.

	This simulator graphs several flight parameters by default, but more are available.  The constants at the top of the FreeRocket class definition control which parameters are to be graphed.  If you encounter a flight parameter which you would like to graph, you can either implement it yourself (the physical processes are not very complicated, and the existing graphs should provide enough example) or submit an issue. 

	All vehicle parameters are modified in the dicts above main(), near the bottom of the file.  Detailed comments for each input are included in the FreeRocket class constructor, so they should appear in your IDE.

Abstract:

	This simulator does not attempt to approximate aerodynamic properties of rockets, only their actual flight characteristics.  This is not intended to be a replacement for OpenRocket or any other full-featured rocket building simulator, only a tool for preliminary analysis.  The simulation possesses 6 degrees of freedom (3-axis translation and 3-axis rotation). 

	As such, you, the user, must supply each property of the vehicle (dry mass, fuel mass, drag coefficient, frontal area, etc.).  There is no attempt to derive any of these values from physical realities, such as the size and density of various structural components or propellants.

	Included in this package are several spreadsheets in OpenDocument format which may help you estimate the mass of various components of your rocket design.

Planned features:

	Thrust vector control and inertial guidance

	Built-in liquid blowdown thrust curve simulation
		This will a re-implementation of the Matlab/Octave script "blowdown.m" into the main Python program
	
Copyright © 2025 John Nguyen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


