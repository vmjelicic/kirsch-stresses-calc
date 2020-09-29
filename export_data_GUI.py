#Librerías generales importadas para GUI
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import filedialog as FileDialog

#Librerias generales del backend
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Importado de scripts de la app
from induced_stresses_calc_backend import BoundaryBox, Geometry_circle, Geometry_ellipse



class ExportDataFrame(ttk.Frame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #Funcionamiento de la pestaña de Exportar Data
        ttk.Label(self, text="Variables to export (use comma and space as separator ', ')").grid(row=0, column=0, sticky = "w")

        self.variable_to_export_entry = ttk.Entry(self, width=30)
        self.variable_to_export_entry.grid(row=1, column=0, pady=5, sticky = (tk.W, tk.E))
        self.variable_to_export_entry.config(state="disabled")

        self.select_variable = tk.IntVar()
        self.select_variable_button = ttk.Checkbutton(self, command=self.active_select_variable, 
                                                      variable=self.select_variable, onvalue=1, offvalue=0)
        self.select_variable_button.grid(row=1, column=3, sticky="e")
        ttk.Label(self, text="Activate selection").grid(row=1, column=4, sticky="w", padx=5)

        ttk.Label(self, text="Folder").grid(row=3, column=0, sticky="w")
        
        self.folder_root_entry = ttk.Entry(self, width=30)
        self.folder_root_entry.grid(row=4, column=0, pady=5, sticky = (tk.W, tk.E))
        self.folder_root_entry.config(state="disabled")

        self.folder_root_button = ttk.Button(self, text="Open", command=self.select_folder_path)
        self.folder_root_button.grid(row=4, column=1, padx=5)

        ttk.Label(self, text="File name").grid(row=5, column=0, sticky="w")

        self.file_name_entry = ttk.Entry(self)
        self.file_name_entry.grid(row=6, column=0, pady=5, sticky=(tk.W, tk.E))

        ttk.Label(self, text=".csv").grid(row=6, column=1, sticky="w", padx=5)

        self.export_button = ttk.Button(self, text="Export data", command=self.export_data) #además cuando termine de exportar debe invocar un pop-up
        self.export_button.grid(row=6, column=5)


    def active_select_variable(self):
        if self.select_variable.get():
            self.variable_to_export_entry.config(state="normal")

        else:
            self.variable_to_export_entry = ttk.Entry(self, width=30)
            self.variable_to_export_entry.grid(row=1, column=0, pady=5, sticky = (tk.W, tk.E))
            self.variable_to_export_entry.config(state="disabled")


    def select_folder_path(self):
        folder_path = FileDialog.askdirectory(title="Select a folder to save file")
        self.folder_root_entry.config(state="normal")
        self.folder_root_entry.delete(0, "end")
        self.folder_root_entry.insert(0, folder_path)
        self.folder_root_entry.config(state="disabled")


    def export_data(self):
        from load_parameters_GUI import boundarybox
        import os

        variables_to_export = self.variable_to_export_entry.get()
        variables_to_export = variables_to_export.split(", ")

        variables_to_export = ["x", "y"] + variables_to_export
        
        folder_path = self.folder_root_entry.get()
        file_name = self.file_name_entry.get()
        file_path = os.path.join(folder_path, file_name + ".csv")

        if self.select_variable.get():
            boundarybox.boundary_df[variables_to_export].to_csv(file_path, sep=",", index=False)
        
        else:
            boundarybox.boundary_df.to_csv(file_path, sep=",", index=False)

        self.variable_to_export_entry.delete(0, "end")
        self.folder_root_entry.config(state="normal")
        self.folder_root_entry.delete(0, "end")
        self.folder_root_entry.config(state="disabled")
        self.file_name_entry.delete(0, "end")

        print("Export successfully")
