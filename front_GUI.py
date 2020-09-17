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


class FrontFrame(ttk.Frame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #Funcionamiento de la pestaña de Inicio
        ttk.Label(self, text="View the instructions on the tabs in order").grid(row=0, column=0, sticky="w", padx=5, pady=10)

        #Creación de notebook
        self.OptionsNotebook = ttk.Notebook(self)

        #Pestaña de Uso de Load Parameters
        self.UseLoadParametersFrame = ttk.Frame(self.OptionsNotebook)
        self.UseLoadParametersFrame.grid(row=0, column=0)

        #Pestaña de Uso de Visualization
        self.UseVisualizationFrame = ttk.Frame(self.OptionsNotebook)
        self.UseVisualizationFrame.grid(row=0, column=0)

        #Pestaña de Uso de Export Data
        self.UseExportDataFrame = ttk.Frame(self.OptionsNotebook)
        self.UseExportDataFrame.grid(row=0, column=0)


        #Información por pestaña
        ########### Pestaña Load Parameters
        self.UseLoadParametersNotebook = ttk.Notebook(self.UseLoadParametersFrame)

        #Pestaña de Uso de Box
        self.box_frame = ttk.Frame(self.UseLoadParametersNotebook)
        self.box_frame.grid(row=0, column=0)
        box_image = tk.PhotoImage(file="./images/box.png")
        box_label = ttk.Label(self.box_frame, image=box_image)
        box_label.image = box_image
        box_label.grid(row=0, column=0)

        #Pestaña de Uso de Geometry_circle
        self.geometry_circle_frame = ttk.Frame(self.UseLoadParametersNotebook)
        self.geometry_circle_frame.grid(row=0, column=0)
        circle_image = tk.PhotoImage(file="./images/circle.png")
        circle_label = ttk.Label(self.geometry_circle_frame, image=circle_image)
        circle_label.image = circle_image
        circle_label.grid(row=0, column=0)

        #Pestaña de Uso de Geometry_ellipse
        self.geometry_ellipse_frame = ttk.Frame(self.UseLoadParametersNotebook)
        self.geometry_ellipse_frame.grid(row=0, column=0)
        ellipse_image = tk.PhotoImage(file="./images/ellipse.png")
        ellipse_label = ttk.Label(self.geometry_ellipse_frame, image=ellipse_image)
        ellipse_label.image = ellipse_image
        ellipse_label.grid(row=0, column=0)

        #Pestaña de Uso de Stresses
        self.stresses_frame = ttk.Frame(self.UseLoadParametersNotebook)
        self.stresses_frame.grid(row=0, column=0)
        stresses_image = tk.PhotoImage(file="./images/hoek-brown.png")
        stresses_label = ttk.Label(self.stresses_frame, image=stresses_image)
        stresses_label.image = stresses_image
        stresses_label.grid(row=0, column=0)

        #Agregado de las pestañas al Notebook
        self.UseLoadParametersNotebook.add(self.box_frame, text="Box", padding=10)
        self.UseLoadParametersNotebook.add(self.geometry_circle_frame, text="Geometry Circle", padding=10)
        self.UseLoadParametersNotebook.add(self.geometry_ellipse_frame, text="Geometry Ellipse", padding=10)
        self.UseLoadParametersNotebook.add(self.stresses_frame, text="Stresses", padding=10)

        #Empaquetado de Notebook - Pestaña Load Parameters
        self.UseLoadParametersNotebook.grid(row=0, column=0)
        ###################################################

        ########### Pestaña Visualization
        # Pestaña Uso de Discretization and Calculate
        visualization_image = tk.PhotoImage(file="./images/visualization.png")
        visualization_label = ttk.Label(self.UseVisualizationFrame, image=visualization_image)
        visualization_label.image = visualization_image
        visualization_label.grid(row=0, column=0)
        ###################################################

        # Pestaña Export Data
        export_data_image = tk.PhotoImage(file="./images/export_data.png")
        export_data_label = ttk.Label(self.UseExportDataFrame, image=export_data_image)
        export_data_label.image = export_data_image
        export_data_label.grid(row=0, column=0)

        #Agregado de las pestañas al Notebook
        self.OptionsNotebook.add(self.UseLoadParametersFrame, text="Load Parameters", padding=10)
        self.OptionsNotebook.add(self.UseVisualizationFrame, text="Visualization", padding=10)
        self.OptionsNotebook.add(self.UseExportDataFrame, text="Export Data", padding=10)

        #Empaquetado notebook
        self.OptionsNotebook.grid(row=1, column=0, pady=5, sticky = "n")



