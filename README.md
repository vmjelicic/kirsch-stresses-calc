# Kirsch's Stresses Calculator

### Scripts list:
- Induced Stresses Calc
- front_GUI
- load_parameters_GUI
- visualization_GUI
- export_data_GUI
- induced_stresses_calc_backend

### Induced Stresses Calc
Is the main script, who use Tkinter for generate a GUI with 4 tabs called: "How to Use" (front_GUI), "Load Parameters" (load_parameters_GUI), "Visualization" (visualization_GUI) and "Export Data" (export_data_GUI). In addition, it has a "dark mode" that modifies the color of the widgets.

### front_GUI
Contain a tab called "How to Use". In this tab, through images (folder "images"), the meaning of the parameters used by the application is shown.

### load_parameters_GUI
Contain a tab called "Load Parameters". In this tab the inputs parameters are loaded, which are divided into: boundary boxes, excavations, stresses fields and strength parameters.

### visualization_GUI
Contain a tab called "Visualization". In this tab the results are shown, where you can select which variable to display from the following: 
- Major principal stress (P)
- Minor principal stress (Q)
- Max shear stress (Tmax)
- Normal stress - X-axis
- Normal stress - Y-axis
- Shear stress- XY-plane
- Strength factor

### export_data_GUI
Contain a tab called "Export Data". In this tab you can export the result of the calculations in a .csv file, where you can select a specific variables or export an entire file.

### induced_stresses_calc_backend
Contain a data treatment and functions used to realize stress calculations, including Kirsch's stresses formulas, rotate functions and overlap stresses functions. Also, this script considers the following objects: "BoundaryBox", "Geometry_circle" and "Geometry_ellipse", where all of calculations are realizing in the BoundaryBox object, while the geometries are used as restrictions for the points generation and as input parameters for each Kirsch's formulae.


##### P.D.: The comments in the scripts are written in Spanish, I must translate them into English.



##### Some views of the app:

<img src="https://github.com/vmyelicich/kirsch-stresses-calc/blob/master/views/Screenshot_6963.png" width="70%" height="70%"/></img>
<img src="https://github.com/vmyelicich/kirsch-stresses-calc/blob/master/views/Screenshot_6964.png" width="70%" height="70%"/></img>
<img src="https://github.com/vmyelicich/kirsch-stresses-calc/blob/master/views/Screenshot_6965.png" width="70%" height="70%"/></img>
<img src="https://github.com/vmyelicich/kirsch-stresses-calc/blob/master/views/Screenshot_6966.png" width="70%" height="70%"/></img>
