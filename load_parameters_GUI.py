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

#Variables globales (no serán importadas, solo están aquí para que no se generen errores de asignación)
boundarybox = 0
geometries = []
P = 0
Q = 0
angle_P = 0
sci = 0
mi = 0
gsi = 0

class LoadParametersFrame(ttk.Frame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #Funcionamiento de la pestaña de Cargar Parámetros
        #Frame Izquierdo
        self.LeftFrame = ttk.Frame(self)
        self.LeftFrame.grid(row=0, column=0, sticky=(tk.N,tk.S), padx=5)
        self.LeftFrame.grid_rowconfigure(2, weight=1)

        #Frame Derecho
        self.RightFrame = ttk.Frame(self)
        self.RightFrame.grid(row=0, column=1, sticky=(tk.N,tk.S))
        self.RightFrame.grid_rowconfigure(2, weight=1)

        #Notebook de opciones para cargado de parámetros
        self.OptionsNotebook = ttk.Notebook(self.RightFrame)

        #Opciones notebook por pestaña

        self.BoxFrame = ttk.Frame(self.OptionsNotebook)
        self.BoxFrame.grid(row=0, column=0, sticky=(tk.N,tk.S))
        self.BoxFrame.grid_rowconfigure(5, weight=1)

        self.dim_x_label = ttk.Label(self.BoxFrame, text="Dim X:", width=10)
        self.dim_x_label.grid(row=0, column=0, sticky="e")
        self.dim_x_entry = ttk.Entry(self.BoxFrame)
        self.dim_x_entry.grid(row=0, column=1)
        self.dim_y_label = ttk.Label(self.BoxFrame, text="Dim Y:", width=10)
        self.dim_y_label.grid(row=1, column=0, sticky="e")
        self.dim_y_entry = ttk.Entry(self.BoxFrame)
        self.dim_y_entry.grid(row=1, column=1)
        self.sep_x_label = ttk.Label(self.BoxFrame, text="Sep X:", width=10)
        self.sep_x_label.grid(row=2, column=0, sticky="e")
        self.sep_x_entry = ttk.Entry(self.BoxFrame)
        self.sep_x_entry.grid(row=2, column=1)
        self.sep_y_label = ttk.Label(self.BoxFrame, text="Sep Y:", width=10)
        self.sep_y_label.grid(row=3, column=0, sticky="e")
        self.sep_y_entry = ttk.Entry(self.BoxFrame)
        self.sep_y_entry.grid(row=3, column=1)
        self.save_box_button = ttk.Button(self.BoxFrame, text="Save", command=self.save_boundarybox)
        self.save_box_button.grid(row=6, column=1)


        self.GeometryCircleFrame = ttk.Frame(self.OptionsNotebook)
        self.GeometryCircleFrame.grid(row=0, column=0, sticky=(tk.N,tk.S))
        self.GeometryCircleFrame.grid_rowconfigure(4, weight=1)

        self.geometry_circle_labels = ["Name:","X Centre:", "Y Centre:", "Radius:"]
        self.geometry_circle_labels_rows = [0, 1, 2, 3]
        for i in range(len(self.geometry_circle_labels)):
            ttk.Label(self.GeometryCircleFrame, text=self.geometry_circle_labels[i], width=10).grid(row=self.geometry_circle_labels_rows[i], column=0, sticky="e")

        self.name_circle_entry = ttk.Entry(self.GeometryCircleFrame)
        self.name_circle_entry.grid(row=0, column=1)

        self.x_centre_circle_entry = ttk.Entry(self.GeometryCircleFrame)
        self.x_centre_circle_entry.grid(row=1, column=1)

        self.y_centre_circle_entry = ttk.Entry(self.GeometryCircleFrame)
        self.y_centre_circle_entry.grid(row=2, column=1)

        self.radius_circle_entry = ttk.Entry(self.GeometryCircleFrame)
        self.radius_circle_entry.grid(row=3, column=1)
        self.save_geometry_circle_button = ttk.Button(self.GeometryCircleFrame, text="Save", command=self.save_geometry_circle)
        self.save_geometry_circle_button.grid(row=6, column=1)



        self.GeometryEllipseFrame = ttk.Frame(self.OptionsNotebook)
        self.GeometryEllipseFrame.grid(row=0, column=0, sticky=(tk.N, tk.S))
        self.GeometryEllipseFrame.grid_rowconfigure(5, weight=1)

        self.geometry_ellipse_labels = ["Name:","X Centre:", "Y Centre:", "Width:", "Height:", "Beta angle:"]
        self.geometry_ellipse_labels_rows = [0, 1, 2, 3, 4, 5]
        for i in range(len(self.geometry_ellipse_labels)):
            ttk.Label(self.GeometryEllipseFrame, text=self.geometry_ellipse_labels[i], width=10).grid(row=self.geometry_ellipse_labels_rows[i], column=0, sticky="e")

        self.name_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.name_ellipse_entry.grid(row=0, column=1)

        self.x_centre_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.x_centre_ellipse_entry.grid(row=1, column=1)

        self.y_centre_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.y_centre_ellipse_entry.grid(row=2, column=1)

        self.width_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.width_ellipse_entry.grid(row=3, column=1)

        self.height_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.height_ellipse_entry.grid(row=4, column=1)

        self.beta_ellipse_entry = ttk.Entry(self.GeometryEllipseFrame)
        self.beta_ellipse_entry.grid(row=5, column=1)
        self.save_geometry_ellipse_button = ttk.Button(self.GeometryEllipseFrame, text="Save", command=self.save_geometry_ellipse)
        self.save_geometry_ellipse_button.grid(row=6, column=1)



        self.StressesFrame = ttk.Frame(self.OptionsNotebook)
        self.StressesFrame.grid(row=0, column=0)

        self.stresses_labels = ["P:", "Q:", "Angle P:", "Sci:", "mi:", "GSI:"]
        self.stresses_labels_rows = [0, 1, 2, 3, 4, 5]
        for i in range(len(self.stresses_labels)):
            ttk.Label(self.StressesFrame, text=self.stresses_labels[i], width=10).grid(row=self.stresses_labels_rows[i], column=0, sticky="e")

        self.P_entry = ttk.Entry(self.StressesFrame)
        self.P_entry.grid(row=0, column=1)

        self.Q_entry = ttk.Entry(self.StressesFrame)
        self.Q_entry.grid(row=1, column=1)

        self.P_angle_entry = ttk.Entry(self.StressesFrame)
        self.P_angle_entry.grid(row=2, column=1)

        self.sci_entry = ttk.Entry(self.StressesFrame)
        self.sci_entry.grid(row=3, column=1)
        self.sci_entry.config(state="disabled")

        self.mi_entry = ttk.Entry(self.StressesFrame)
        self.mi_entry.grid(row=4, column=1)
        self.mi_entry.config(state="disabled")

        self.gsi_entry = ttk.Entry(self.StressesFrame)
        self.gsi_entry.grid(row=5, column=1)
        self.gsi_entry.config(state="disabled")
        self.save_stresses_button = ttk.Button(self.StressesFrame, text="Save", command=self.save_stresses)
        self.save_stresses_button.grid(row=6, column=1)

        self.strength_option = tk.IntVar()

        self.strength_button = ttk.Checkbutton(self.StressesFrame, command=self.active_strength_inputs, 
                                                variable=self.strength_option, onvalue=1, offvalue=0)
        self.strength_button.grid(row=3, column=2)
        self.strength_label = ttk.Label(self.StressesFrame, text="Strength parameters")
        self.strength_label.grid(row=3, column=3)

        
        #Widgets frame izquierdo
        self.saved_parameters_label = ttk.Label(self.LeftFrame, text="History")
        self.saved_parameters_label.grid(row=0, column=0, sticky="w", pady=5)

        self.saved_parameters_text = tk.Text(self.LeftFrame)
        self.saved_parameters_text.grid(row=1, column=0, pady=5, sticky=(tk.N,tk.S))
        self.saved_parameters_text.config(width=70, height=17, padx=10, pady=5, state="disabled")

        #Widgets frame derecho
        self.load_parameters_label = ttk.Label(self.RightFrame, text="Load parameters")
        self.load_parameters_label.grid(row=0, column=0, sticky="w", pady=5)

        #Opciones Notebook
        self.OptionsNotebook.add(self.BoxFrame, text="Box", padding=10)
        self.OptionsNotebook.add(self.GeometryCircleFrame, text="Geometry Circle", padding=10)
        self.OptionsNotebook.add(self.GeometryEllipseFrame, text="Geometry Ellipse", padding=10)
        self.OptionsNotebook.add(self.StressesFrame, text="Stresses", padding=10)

        #Modificaciones de historial de frame izquierdo
        self.delete_parameter_button = ttk.Button(self.RightFrame, text="Delete Parameter", width=15) #command=
        self.delete_parameter_button.grid(row=2, column=0, sticky="s", pady=5)
        self.clear_history_button = ttk.Button(self.RightFrame, text="Clear History", command=self.clear_history, width=15)
        self.clear_history_button.grid(row=3, column=0, sticky="s", pady=5)

        #Empaquetado notebook frame derecho
        self.OptionsNotebook.grid(row=1, column=0, pady=5)

    def active_strength_inputs(self):
        if (self.strength_option.get()):
            self.sci_entry.config(state="normal")
            self.mi_entry.config(state="normal")
            self.gsi_entry.config(state="normal")
        else:
            self.sci_entry = ttk.Entry(self.StressesFrame)
            self.sci_entry.grid(row=3, column=1)
            self.sci_entry.config(state="disabled")
            self.mi_entry = ttk.Entry(self.StressesFrame)
            self.mi_entry.grid(row=4, column=1)
            self.mi_entry.config(state="disabled")
            self.gsi_entry = ttk.Entry(self.StressesFrame)
            self.gsi_entry.grid(row=5, column=1)
            self.gsi_entry.config(state="disabled")


    #Funciones de los distintos botones
    
    def save_boundarybox(self):
        dim_x = float(self.dim_x_entry.get())
        dim_y = float(self.dim_y_entry.get())
        sep_x = float(self.sep_x_entry.get())
        sep_y = float(self.sep_y_entry.get())

        #Creación de objeto BoundaryBox y guardado
        global boundarybox
        boundarybox = BoundaryBox(dim_x=dim_x, dim_y=dim_y, sep_x=sep_x, sep_y=sep_y)
        

        #Escritura de registro en el widget text
        self.saved_parameters_text.config(state="normal")
        self.saved_parameters_text.insert(tk.END, "BoundaryBox(dim_x={}, dim_y={}, sep_x={}, sep_y={})\n".format(dim_x, dim_y, sep_x, sep_y))
        self.saved_parameters_text.config(state="disabled")

        #Limpieza de entries
        self.dim_x_entry.delete(0, "end")
        self.dim_y_entry.delete(0, "end")
        self.sep_x_entry.delete(0, "end")
        self.sep_y_entry.delete(0, "end")


    def save_geometry_circle(self):
        name = self.name_circle_entry.get()
        x_centre = float(self.x_centre_circle_entry.get())
        y_centre = float(self.y_centre_circle_entry.get())
        radius = float(self.radius_circle_entry.get())

        #Creación de objeto Geometry_circle
        new_geometry_circle = Geometry_circle(name=name, x_centre=x_centre, y_centre=y_centre, radius=radius)
        #Guardado de objeto
        global geometries
        geometries.append(new_geometry_circle)

        #Escritura de registro en el widget text
        self.saved_parameters_text.config(state="normal")
        self.saved_parameters_text.insert(tk.END, "Geometry_circle(name={}, x_centre={}, y_centre={}, radius={})\n".format(name, x_centre, y_centre, radius))
        self.saved_parameters_text.config(state="disabled")

        #Limpieza de entries
        self.name_circle_entry.delete(0, "end")
        self.x_centre_circle_entry.delete(0, "end")
        self.y_centre_circle_entry.delete(0, "end")
        self.radius_circle_entry.delete(0, "end")


    def save_geometry_ellipse(self):
        name = self.name_ellipse_entry.get()
        x_centre = float(self.x_centre_ellipse_entry.get())
        y_centre = float(self.y_centre_ellipse_entry.get())
        width = float(self.width_ellipse_entry.get())
        height = float(self.height_ellipse_entry.get())
        beta = float(self.beta_ellipse_entry.get())

        #Creación del objeto Geometry_ellipse
        new_geometry_ellipse = Geometry_ellipse(name=name, x_centre=x_centre, y_centre=y_centre, width=width, height=height, beta=beta)
        #Guardado de objeto
        global geometries
        geometries.append(new_geometry_ellipse)

        #Escritura de registro en el widget text
        self.saved_parameters_text.config(state="normal")
        self.saved_parameters_text.insert(tk.END, "Geometry_ellipse(name={}, x_centre={}, y_centre={}, width={}, height={}, beta={})\n".format(name, x_centre, y_centre, width, height, beta))
        self.saved_parameters_text.config(state="disabled")

        #Limpieza de entries
        self.name_ellipse_entry.delete(0, "end")
        self.x_centre_ellipse_entry.delete(0, "end")
        self.y_centre_ellipse_entry.delete(0, "end")
        self.width_ellipse_entry.delete(0, "end")
        self.height_ellipse_entry.delete(0, "end")
        self.beta_ellipse_entry.delete(0, "end")


    def save_stresses(self):
        #Indicador de parametros de resistencia
        option = self.strength_option.get()

        #Guardado de valores de esfuerzos
        global P, Q, angle_P
        P = float(self.P_entry.get()) ##### Se guarda en las variables globales
        Q = float(self.Q_entry.get())
        angle_P = float(self.P_angle_entry.get())

        #Escritura de registro en el widget text
        self.saved_parameters_text.config(state="normal")
        self.saved_parameters_text.insert(tk.END, "P={}, Q={}, angle_P={}\n".format(P, Q, angle_P))
        self.saved_parameters_text.config(state="disabled")

        #Limpieza de entries
        self.P_entry.delete(0, "end")
        self.Q_entry.delete(0, "end")
        self.P_angle_entry.delete(0, "end")

        #["P:", "Q:", "Angle P:", "Sci:", "mi:", "GSI:"]
        #Si se incorporarán parámetros de resistencia
        if option == True:
            global sci, mi, gsi
            #Guardado de valores de resistencia
            sci = float(self.sci_entry.get()) ##### Se guarda en las variables globales
            mi = float(self.mi_entry.get())
            gsi = float(self.gsi_entry.get())

            #Escritura de registro en el widget text
            self.saved_parameters_text.config(state="normal")
            self.saved_parameters_text.insert(tk.END, "Sci={}, mi={}, GSI={}\n".format(sci, mi, gsi))
            self.saved_parameters_text.config(state="disabled")

            #Limpieza de entries
            self.sci_entry.delete(0, "end")
            self.mi_entry.delete(0, "end")
            self.gsi_entry.delete(0, "end")


    def clear_history(self):
        global boundarybox, geometries

        #Eliminación de boundaries
        boundarybox = 0

        #Eliminación de geometrias
        geometries = []

        #Reset de esfuerzos/resistencia
        global P, Q, angle_P, sci, mi, gsi
        P = 0
        Q = 0
        angle_P = 0
        sci = 0
        mi = 0
        gsi = 0

        #Reset de registro del widget text
        self.saved_parameters_text.config(state="normal")
        self.saved_parameters_text.delete(1.0, "end")
        self.saved_parameters_text.config(state="disabled")

        print("Parameters deleted")