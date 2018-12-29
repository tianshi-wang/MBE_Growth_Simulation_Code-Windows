# Source Files of the MBE Growth Simulation Code (Windows)

# Table of content
1. [Description](README.md#Description)
2. [Approach](README.md#approach)
3. [User Interface](REAME.md#user-interface)
4. [Code Structure](README.md#code-structure)
5. [Contacts](README.md#contacts)

# Description
This is the source file package of the MBE Growth Simulation Code for Windows users. For Windows users, please find the source files and compiled codes using this [link](https://github.com/tianshi-wang/MBE_Growth_Simulation_Code-Mac).

The code can simulate growth morphology at different conditions for molecular beam epitaxy (MBE). The calculation is based on kinetic Monte Carlo method. The interactive 3D rending is excellent for education purpose.

# Approach
The code is written using C++ with OpenGL and QT5. To revise or recompile the code, first download [QT5](https://www1.qt.io/download-open-source/?hsCtaTracking=f977210e-de67-475f-a32b-65cec207fd03%7Cd62710cd-e1db-46aa-8d4d-2f1c1ffdacea#section-2). Then, open the downloaded source files (**MBEcodeSourceFiles_MAC.zip**) with QT5. Make sure your framework has the OpenGL and GLUT libraries. 

# User Interface 

Code layout is shown below. 

<image height=400 src="https://github.com/tianshi-wang/MBE_Growth_Simulation_Code-Mac/blob/master/screenshots/Initial_state.png">
<br>

By adjusting parameters, the following growth morphologies can be observed.

<image height=220 src="https://github.com/tianshi-wang/MBE_Growth_Simulation_Code-Mac/blob/master/screenshots/morphologies.png">
<br>


# Code Structure
- Icon.icon: Code icon. 
- MBEcodeProject.pro and MBEcodeProject.pro.user: QT5 profile including required libraries
- griddata.cpp and griddata.h: Setting different grid dimensions and related variables. 
- main.cpp: Main code.
- myglwidget.cpp, myglwidget.h: OpenGL 3D rending
- results.cpp, results.h: Output Monte-Carlo calculatin results
- solver.cpp, solver.h: KMC solver with thread control
- window.cpp, window.h, and window.ui: QT5 window and widges

# Contacts
For any question, please contact us: Tianshi Wang (tswang@udel.edu) & Wei Li (Verali@udel.edu).
