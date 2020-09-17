#Librerías generales importadas para GUI
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

#Librerias generales del backend
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Importado de scripts de la app
from induced_stresses_calc_backend import BoundaryBox, Geometry_circle, Geometry_ellipse

#Funciones auxiliares para rotación
def transform_xt(x, y, theta):
    return x*np.cos(theta) - y*np.sin(theta)

def transform_yt(x, y, theta):
    return x*np.sin(theta) + y*np.cos(theta)



class VisualizationFrame(ttk.Frame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #Funcionamiento de la pestaña de Visualización

        self.PlotFrame = ttk.Frame(self)
        self.PlotFrame.grid(row=0, column=0)

        self.OptionsFrame = ttk.Frame(self)
        self.OptionsFrame.grid(row=0, column=1, sticky=(tk.N,tk.S), padx=10, pady=10)
        self.OptionsFrame.grid_rowconfigure(9, weight=1)

        #Panel de opciones
        self.discretization_button = ttk.Button(self.OptionsFrame, text="Discretization", command=self.discretization, width=12)
        self.discretization_button.grid(row=0, column=0, sticky="n", pady=5)
        
        self.calculate_button = ttk.Button(self.OptionsFrame, text="Calculate", command=self.calculate, width=12)
        self.calculate_button.grid(row=1, column=0, sticky="n", pady=5)

        ttk.Label(self.OptionsFrame, text="Select variable to show").grid(row=2, column=0, sticky="s")

        self.option_list = ["Variable", "P", "Q", "Tmax", "Sx", "Sy", "Txy", "Strength_factor"]
        self.selected_variable = tk.StringVar()
        self.variable_optionmenu = ttk.OptionMenu(self.OptionsFrame, self.selected_variable, *self.option_list)
        self.variable_optionmenu.grid(row=3, column=0, sticky="s")

        ttk.Label(self.OptionsFrame, text="Activate plot options").grid(row=4, column=0)
        self.plot_options = tk.IntVar()
        self.plot_options_button = ttk.Checkbutton(self.OptionsFrame, command=self.active_plot_options, 
                                                    variable=self.plot_options, onvalue=1, offvalue=0)
        self.plot_options_button.grid(row=4, column=1, sticky="w")
        
        ttk.Label(self.OptionsFrame, text="Min value:").grid(row=5, column=0)
        self.min_value_plot_entry = ttk.Entry(self.OptionsFrame)
        self.min_value_plot_entry.grid(row=5, column=1)
        self.min_value_plot_entry.config(state="disabled")

        ttk.Label(self.OptionsFrame, text="Max value:").grid(row=6, column=0)
        self.max_value_plot_entry = ttk.Entry(self.OptionsFrame)
        self.max_value_plot_entry.grid(row=6, column=1)
        self.max_value_plot_entry.config(state="disabled")

        ttk.Label(self.OptionsFrame, text="Point size:").grid(row=7, column=0)
        self.point_size_entry = ttk.Entry(self.OptionsFrame)
        self.point_size_entry.grid(row=7, column=1)
        self.point_size_entry.config(state="disabled")

        self.select_plot_button = ttk.Button(self.OptionsFrame, text="Generate Plot", command=self.select_plot, width=12)
        self.select_plot_button.grid(row=8, column=0, sticky="s", pady=10)

        self.progress_label = ttk.Label(self.OptionsFrame, text="Stand-by", width=26)
        self.progress_label.grid(row=9, column=0, sticky="se")

        
        self.progress_bar = ttk.Progressbar(self.OptionsFrame, orient=tk.HORIZONTAL, length=100, mode="determinate")
        self.progress_bar.grid(row=9, column=1, sticky="s")

        self.update_bar(value=0)

        self.graph = self.matplot_canvas(x_plot=0, y_plot=0)
        self.toolbar = NavigationToolbar2Tk(self.graph, self.PlotFrame)
        self.toolbar.update()
        self.graph.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        

    #Plot de gráficos
    def matplot_canvas(self, x_plot, y_plot, c_plot=None, point_size=1, label_plot=None, min_value=None, max_value=None, geometries=None):
        f = Figure(figsize=(7,5), dpi=100)
        a = f.add_subplot(111)
        graph = a.scatter(x=x_plot, y=y_plot, c=c_plot, cmap="coolwarm", s=point_size)

        if geometries is not None:
            for geometry in geometries:
                if geometry.__class__.__name__ == "Geometry_circle":
                    x_geom_upper = [x + geometry.x_centre - geometry.radius for x in np.arange(0, 2*geometry.radius+1,1)]
                    x_geom_lower = [x + geometry.x_centre - geometry.radius for x in np.arange(0, 2*geometry.radius+1,1)]
                    y_geom_upper_contour = [np.sqrt(geometry.radius**2 - (x - geometry.x_centre)**2)+geometry.y_centre for x in x_geom_upper]
                    y_geom_lower_contour = [-np.sqrt(geometry.radius**2 - (x - geometry.x_centre)**2)+geometry.y_centre for x in x_geom_lower]

                elif geometry.__class__.__name__ == "Geometry_ellipse":
                    #Generar puntos en el sistema de referencia (x,y)
                    x_geom = [x + geometry.x_centre - (geometry.W/2) for x in np.arange(0, geometry.W + 1, 1)]
                    y_geom_u = [(geometry.H/2)*np.sqrt(1-((x-geometry.x_centre)/(geometry.W/2))**2) + geometry.y_centre for x in x_geom]
                    y_geom_l = [-(geometry.H/2)*np.sqrt(1-((x-geometry.x_centre)/(geometry.W/2))**2) + geometry.y_centre for x in x_geom]
                
                    x_geom_upper = [1 for x in range(len(x_geom))]
                    y_geom_upper_contour = [1 for x in range(len(x_geom))]
                    x_geom_lower = [1 for x in range(len(x_geom))]
                    y_geom_lower_contour = [1 for x in range(len(x_geom))]
                    #Rotar puntos a sistema de referencia (x',y')
                    for i in range(len(x_geom)):
                        x_geom_upper[i] = transform_xt(x = x_geom[i]-geometry.x_centre, y = y_geom_u[i]-geometry.y_centre, theta = geometry.beta) + geometry.x_centre
                        y_geom_upper_contour[i] = transform_yt(x = x_geom[i]-geometry.x_centre, y = y_geom_u[i]-geometry.y_centre, theta = geometry.beta) + geometry.y_centre
                        x_geom_lower[i] = transform_xt(x = x_geom[i]-geometry.x_centre, y = y_geom_l[i]-geometry.y_centre, theta = geometry.beta) + geometry.x_centre
                        y_geom_lower_contour[i] = transform_yt(x = x_geom[i]-geometry.x_centre, y = y_geom_l[i]-geometry.y_centre, theta = geometry.beta) + geometry.y_centre


                a.plot(x_geom_upper, y_geom_upper_contour, "k")
                a.plot(x_geom_lower, y_geom_lower_contour, "k")

        if (min_value is not None) and (max_value is not None):
            f.colorbar(graph, label=label_plot, pad=0.05).mappable.set_clim(min_value, max_value)

        else:
            f.colorbar(graph, label=label_plot, pad=0.05)

        a.grid(color='white', linestyle="--", linewidth=0.5)

        canvas = FigureCanvasTkAgg(f, self.PlotFrame)
        canvas.draw()
        
        return canvas

    #Actualización de barra de progreso
    #Adaptar esta barra a los cálculos
    def update_bar(self, value):
        self.progress_bar["value"] = value
        self.OptionsFrame.update_idletasks()

    #Funciones de los botones
    def discretization(self):
        #Variables globales importadas de load_parameters
        from load_parameters_GUI import boundarybox, geometries, P, Q, angle_P, sci, mi, gsi
        #Creación de puntos
        x_point, y_point = boundarybox.mesh_points(geometries)
        
        self.graph.get_tk_widget().pack_forget()
        self.toolbar.destroy()

        self.graph = self.matplot_canvas(x_plot=x_point, y_plot=y_point, geometries=geometries)
        self.toolbar = NavigationToolbar2Tk(self.graph, self.PlotFrame)
        self.toolbar.update()
        self.graph.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)


    def calculate(self):
        #Variables globales importadas de load_parameters
        from load_parameters_GUI import boundarybox, geometries, P, Q, angle_P, sci, mi, gsi

        #Realización de cálculos de esfuerzos, actualización de barra de progreso y label
        #Opciones por si se incluyen calculos de resistencia o no
        strength_option = (sci > 0) and (mi > 0) and (100 >= gsi > 0)

        if strength_option is True:
            self.progress_label.config(text="Calculating Kirsch's stresses")
            self.update_bar(value=0)

            boundarybox.calculate_kirsch_stresses(P=P, Q=Q, angle_P=angle_P)
            self.progress_label.config(text="Calculating Sx, Sy and Txy")
            self.update_bar(value=25)

            boundarybox.calculate_stress_overlap()
            self.progress_label.config(text="Calculating principal stresses")
            self.update_bar(value=50)

            boundarybox.calculate_principal_stresses()
            self.progress_label.config(text="Calculating strength factor")
            self.update_bar(value=75)

            boundarybox.calculate_strength_factor(sci=sci, mi=mi, gsi=gsi)
            self.update_bar(value=100)
            self.progress_label.config(text="Calculation complete")

        #Realización de cálculos de resistencia
        
        else:
            self.progress_label.config(text="Calculating Kirsch's stresses")
            self.update_bar(value=0)

            boundarybox.calculate_kirsch_stresses(P=P, Q=Q, angle_P=angle_P)
            self.progress_label.config(text="Calculating Sx, Sy and Txy")
            self.update_bar(value=33)

            boundarybox.calculate_stress_overlap()
            self.progress_label.config(text="Calculating principal stresses")
            self.update_bar(value=66)

            boundarybox.calculate_principal_stresses()
            self.update_bar(value=100)
            self.progress_label.config(text="Calculation complete")
        #Hasta acá estarían los calculos listos y sus resultados en el boundary.boundary_df

    def active_plot_options(self):
        if (self.plot_options.get()):
            self.min_value_plot_entry.config(state="normal")
            self.max_value_plot_entry.config(state="normal")
            self.point_size_entry.config(state="normal")

        else:
            self.min_value_plot_entry = ttk.Entry(self.OptionsFrame)
            self.min_value_plot_entry.grid(row=5, column=1)
            self.min_value_plot_entry.config(state="disabled")
            self.max_value_plot_entry = ttk.Entry(self.OptionsFrame)
            self.max_value_plot_entry.grid(row=6, column=1)
            self.max_value_plot_entry.config(state="disabled")
            self.point_size_entry = ttk.Entry(self.OptionsFrame)
            self.point_size_entry.grid(row=7, column=1)
            self.point_size_entry.config(state="disabled")


    def select_plot(self):
        from load_parameters_GUI import boundarybox, geometries
        x_point = boundarybox.boundary_df["x"]
        y_point = boundarybox.boundary_df["y"]
        option_plot = self.selected_variable.get()
        c_point = boundarybox.boundary_df[option_plot]
        
        self.graph.get_tk_widget().pack_forget()
        self.toolbar.destroy()

        ######Revisar estas opciones, no funciona ######### -> Revisar con la función matplot_canvas también
        if (self.plot_options.get()):
            min_value = float(self.min_value_plot_entry.get())
            max_value = float(self.max_value_plot_entry.get())
            point_size = float(self.point_size_entry.get())

        else:
            min_value = None
            max_value = None
            point_size = 1

        self.graph = self.matplot_canvas(x_plot=x_point, y_plot=y_point, 
                                        geometries=geometries, c_plot=c_point, 
                                        label_plot=option_plot, min_value=min_value, 
                                        max_value=max_value, point_size=point_size)

        self.toolbar = NavigationToolbar2Tk(self.graph, self.PlotFrame)
        self.toolbar.update()
        self.graph.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        #self.option_list = ["Variable", "P", "Q", "Tmax", "Sx", "Sy", "Txy"]
        #self.selected_variable = tk.StringVar()
        
        #Usar funciones del backend, pero que junto con actualizar el dataframe vayan retornando valores para los plots
        #Puedo acceder al dataframe boundarybox.boundary_df["variable"]
        #Tambien esto debe coincidir con la variable OptionMenu, para mostrar cada plot

