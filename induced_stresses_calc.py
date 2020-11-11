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
from front_GUI import FrontFrame
from load_parameters_GUI import LoadParametersFrame
from visualization_GUI import VisualizationFrame
from export_data_GUI import ExportDataFrame
from induced_stresses_calc_backend import BoundaryBox, Geometry_circle, Geometry_ellipse

#Variables globales importadas de load_parameters
from load_parameters_GUI import boundarybox, geometries, P, Q, angle_P, sci, mi, gsi


class Application(ttk.Frame):
    
    def __init__(self, main_window):
        super().__init__(main_window)
        main_window.title("Induced Stresses Calc")
        main_window.iconbitmap("./images/icono.ico")
        main_window.resizable(1,1)

        # Crear el estilo de los widgets
        self.style_notebook = ttk.Style()
        self.style_notebook.configure("TNotebook", foreground="black", background="gray94")
        self.style_frames = ttk.Style()
        self.style_frames.configure("TFrame", foreground="black", background="gray94")
        self.style_labels = ttk.Style()
        self.style_labels.configure("TLabel", foreground="black", background="gray94")
        self.style_entries = ttk.Style()
        self.style_entries.configure("TEntry", foreground="black", background="gray94")
        self.style_checkbuttons = ttk.Style()
        self.style_checkbuttons.configure("TCheckbutton", background="gray94")
        self.style_text = ttk.Style()
        self.style_text.configure("TText", foreground="black", background="gray94")
        self.style_progressbar = ttk.Style()


        # Crear el panel de pestañas.
        self.notebook = ttk.Notebook(self, style="TNotebook")

        # Crear el contenido de cada una de las pestañas.
        self.front_label = FrontFrame(self.notebook)
        self.load_parameters_label = LoadParametersFrame(self.notebook)
        self.visualization_label = VisualizationFrame(self.notebook)
        self.export_data_label = ExportDataFrame(self.notebook)


        # Añadirlas al panel con su respectivo texto.
        self.notebook.add(self.front_label, text="How to Use", padding=5)
        self.notebook.add(self.load_parameters_label, text="Load Parameters", padding=5)
        self.notebook.add(self.visualization_label, text="Visualization", padding=5)
        self.notebook.add(self.export_data_label, text="Export Data", padding=5)
        self.notebook.grid(row=0, column=0, padx=10, pady=10)
        

        self.pack()

        self.dark_mode_var = tk.IntVar()
        ttk.Label(self, text="Dark mode").grid(row=0, column=1, sticky="ne")
        ttk.Checkbutton(self, command=self.dark_mode, variable=self.dark_mode_var, onvalue=1, offvalue=0).grid(row=0, column=2, sticky="ne")


    def dark_mode(self):
        if (self.dark_mode_var.get()):
            self.style_notebook.configure("TNotebook", foreground="white", background="gray20")
            self.style_frames.configure("TFrame", foreground="white", background="gray20")
            self.style_labels.configure("TLabel", foreground="white", background="gray20")    
            self.style_entries.configure("TEntry", foreground="black", background="gray20")
            self.style_checkbuttons.configure("TCheckbutton", background="gray20")
            self.style_text.configure("TText", foreground="white", background="gray20")    

        else:
            self.style_notebook.configure("TNotebook", foreground="black", background="gray94")
            self.style_frames.configure("TFrame", foreground="black", background="gray94")
            self.style_labels.configure("TLabel", foreground="black", background="gray94")    
            self.style_entries.configure("TEntry", foreground="black", background="gray94")
            self.style_checkbuttons.configure("TCheckbutton", background="gray94")
            self.style_text.configure("TText", foreground="black", background="gray94") 


#Ejecución de estas lineas solo si este script es el principal (no importado)
if __name__ == "__main__":
    main_window = tk.Tk()
    app = Application(main_window)
    app.mainloop()
