# Librerías Importadas
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Definición de clases
# Geometría circular
class Geometry_circle:
    def __init__(self, name, x_centre, y_centre, radius):
        self.name = name
        self.x_centre = x_centre
        self.y_centre = y_centre
        self.radius = radius
        print("Geometry (circle) parameters saved")

# Geometría Elíptica
class Geometry_ellipse:
    def __init__(self, name, x_centre, y_centre, W, H, beta):
        #Comprobación de parámetros correctos
        if W == H:
            print("Invalid ellipse dimensions")
        else:
            #Id de geometría
            self.name = name
            #Coordenadas del centro
            self.x_centre = x_centre
            self.y_centre = y_centre
            #Extensión de ejes de elipse 
            #Semiejes son W/2 y H/2
            self.W = W
            self.H = H
            #Ángulo de inclinación de la elipse
            self.beta = np.radians(beta)
            print("Geometry (ellipse) parameters saved")

# Boundary box, o caja en la cual se realizarán los cálculos
class BoundaryBox:
    def __init__(self, dim_x, dim_y, sep_x, sep_y):
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.sep_x = sep_x
        self.sep_y = sep_y
        self.boundary_df = pd.DataFrame({"x":[],"y":[]})
        print("Boundary box parameters saved")

    # Rotación de sistema de referencia de (x,y) a (x',y') en función de ángulo theta (positivo en sentido antihorario)
    def transform_x(self, x, y, theta):
        return x*np.cos(theta) + y*np.sin(theta)

    def transform_y(self, x, y, theta):
        return -x*np.sin(theta) + y*np.cos(theta)

    #Rotación de sistema de referencia de (x',y') a (x,y) en función de ángulo theta (positivo en sentido horario)
    def transform_xt(self, x, y, theta):
        return x*np.cos(theta) - y*np.sin(theta)

    def transform_yt(self, x, y, theta):
        return x*np.sin(theta) + y*np.cos(theta)

    # Método de discretización o generación de puntos para realizar cálculos
    def mesh_points(self, geometries):
        # Guardado de geometrías
        self.geometries = geometries

        # Listas para guardar resultados
        x = []
        y = []
        # Contador de cumplimiento de restricción por geometría
        out_of_geometry = True

        # Saltos definidos en X e Y (sep_x, sep_y)
        for i in np.arange(-self.dim_x, self.dim_x, self.sep_x):
            for j in np.arange(-self.dim_y, self.dim_y, self.sep_y):
                
                # Condición para que no se generen puntos dentro de las geometrías, para cada geometría
                for geometry in self.geometries:
                    if geometry.__class__.__name__ == "Geometry_circle":
                        if ((i-geometry.x_centre)**2 + (j-geometry.y_centre)**2) < geometry.radius**2:
                            out_of_geometry = False

                    elif geometry.__class__.__name__ == "Geometry_ellipse":
                        #Transformación de sistema de coordenadas (o rotación) para coordenadas (x-h, y-k)
                        x_tf = self.transform_x(x = i-geometry.x_centre, y = j-geometry.y_centre, theta = geometry.beta) + geometry.x_centre
                        y_tf = self.transform_y(x = i-geometry.x_centre, y = j-geometry.y_centre, theta = geometry.beta) + geometry.y_centre

                        if (((x_tf-geometry.x_centre)**2)/((geometry.W/2)**2) + ((y_tf-geometry.y_centre)**2)/((geometry.H/2)**2)) < 1:
                            out_of_geometry = False

                if out_of_geometry:
                    x.append(i)
                    y.append(j)

                # Reseteo de contador al finalizar la asignación
                out_of_geometry = True 

        # Guardado de información en el dataframe
        self.boundary_df["x"] = x
        self.boundary_df["y"] = y
        print("Meshing performed")

        return x, y

    # Métodos para cálculos de esfuerzos
    # Cálculo de "r" por punto
    def kirsch_r(self, x_point, y_point, x_centre_geometry, y_centre_geometry):
        return np.sqrt((x_point - x_centre_geometry)**2 + (y_point - y_centre_geometry)**2)

    # Cálculo de ángulo "theta" por punto (en radianes)
    def kirsch_theta(self, x_point, y_point, x_centre_geometry, y_centre_geometry):
        # Evitar divisiones por cero (cuando theta = 90°)
        if (x_centre_geometry - x_point) == 0:
            return np.pi/2
        else:
            return np.arctan((y_point - y_centre_geometry) / (x_point - x_centre_geometry))

    # Funciones de cálculo de esfuerzos inducidos (circunferencia)
    def kirsch_sigma_r(self, sx, sy, txy, r, a, theta):
        return ((sx+sy)/2) * (1-(a/r)**2) + (1+3*(a/r)**4-4*(a/r)**2) * (((sx-sy)/2)*np.cos(2*theta) + txy*np.sin(2*theta))

    def kirsch_sigma_theta(self, sx, sy, txy, r, a, theta):
        return ((sx+sy)/2) * (1+(a/r)**2) - (1+3*(a/r)**4) * (((sx-sy)/2)*np.cos(2*theta) + txy*np.sin(2*theta))

    def kirsch_tau_r_theta(self, sx, sy, txy, r, a, theta):
        return (1-3*(a/r)**4+2*(a/r)**2) * (((sy-sx)/2)*np.sin(2*theta) + txy*np.cos(2*theta))

    # Esfuerzos inducidos (Elipse)
        # Cálculo de variable auxiliar "e0"
    def kirsch_ellipse_e0(self, W, H):
        return (W+H)/(W-H)

    # Cálculo de variable auxiliar "b"
    def kirsch_ellipse_b(self, W, H, x1, z1):
        return (4*(x1**2 + z1**2))/(W**2 - H**2)

    # Cálculo de variable auxiliar "d"
    def kirsch_ellipse_d(self, W, H, x1, z1):
        return (8*(x1**2 - z1**2))/(W**2 - H**2) - 1

    # Cálculo de variable auxiliar "u"
    def kirsch_ellipse_u(self, b, e0, d):
        return b + (e0/abs(e0)) * (b**2 - d)**0.5

    # Cálculo de variable auxiliar "e"
    def kirsch_ellipse_e(self, e0, u):
        return u + (e0/abs(e0)) * (u**2 - 1)**0.5

    # Cálculo de variable auxiliar "phi"
    def kirsch_ellipse_phi(self, e, x1, z1):
        if ((e-1) == 0) or (x1 == 0):
            return np.pi/2
        else:
            return np.arctan(((e+1)/(e-1))*(z1/x1))

    # Cálculo de variable auxiliar "theta"
    def kirsch_ellipse_theta(self, e, x1, z1):
        if ((e-1) == 0) or (x1 == 0):
            return np.pi/2
        else:
            return np.arctan((((e+1)/(e-1))**2)*(z1/x1))

    # Cálculo de variable auxiliar "C"
    def kirsch_ellipse_C(self, e0, e):
        return 1 - e*e0

    # Cálculo de variable auxiliar "J"
    def kirsch_ellipse_J(self, e, phi):
        return 1 + e**2 - 2*e*np.cos(2*phi)

    # Cálculo de esfuerzo radial
    def kirsch_ellipse_sigma_r(self, p, e0, e, C, J, K, phi, beta):
        return ((p*(e0-e))/(J**2)) * ((1+K)*(e**2-1)*(C/(2*e0))+(1-K)*(((J/2)*(e-e0)+C*e)*np.cos(2*(phi+beta))-C*np.cos(2*beta)))

    # Cálculo de esfuerzo tangencial
    def kirsch_ellipse_sigma_theta(self, p, e0, e, J, K, phi, beta, sigma_r):
        return (p/J) * ((1+K)*(e**2-1)+2*(1-K)*e0*(e*np.cos(2*(phi+beta))-np.cos(2*beta))) - sigma_r

    # Cálculo de esfuerzo cortante
    def kirsch_ellipse_tau_r_theta(self, p, e0, e, C, J, K, phi, beta):
        return ((p*(e0-e))/(J**2)) * ((1+K)*((C*e)/e0)*np.sin(2*phi)+(1-K)*(e*(e0+e)*np.sin(2*beta)+e*np.sin(2*(phi-beta))-((J/2)*(e0+e)+e**2*e0)*np.sin(2*(phi+beta))))

    # Funciones de rotación de esfuerzos 2d
    # Cálculo de Sx
    def rotate_sigma_x(self, sx, sy, txy, theta):
        return ((sx+sy)/2) + ((sx-sy)/2) * np.cos(2*theta) + txy*np.sin(2*theta)

    # Cálculo de Sy
    def rotate_sigma_y(self, sx, sy, txy, theta):
        return ((sx+sy)/2) - ((sx-sy)/2) * np.cos(2*theta) - txy*np.sin(2*theta)

    # Cálculo de Txy
    def rotate_tau_xy(self, sx, sy, txy, theta):
        return txy*np.cos(2*theta) - ((sx-sy)/2) * np.sin(2*theta)

    # Funciones de esfuerzos principales
    # Cálculo de P
    def principal_stress_P(self, sx, sy, txy):
        return ((sx+sy)/2) + np.sqrt(((sx-sy)/2)**2 + txy**2)

    # Cálculo de Q
    def principal_stress_Q(self, sx, sy, txy):
        return ((sx+sy)/2) - np.sqrt(((sx-sy)/2)**2 + txy**2)

    # Cálculo de Tmax
    def principal_stress_max_shear(self, sx, sy, txy):
        return np.sqrt(((sx-sy)/2)**2 + txy**2)

    # Funciones de cálculo de resistencia
    # Cálculo de envolvente de falla de Hoek - Brown
    def hoek_brown_S1(self, s3, sci_rock, mi_rock, gsi_rockmass):
        mb = mi_rock * np.exp((gsi_rockmass-100)/28)
        s = np.exp((gsi_rockmass-100)/9)
        a = 0.5 + (1/6)*(np.exp(-gsi_rockmass/15)-np.exp(-20/3))
        if s3 < 0:
            return 0
        else:
            return s3 + sci_rock*((mb/sci_rock)*s3 + s)**a

    # Cálculo de factor de seguridad
    def strength_factor(self, s1_induced, s1_rockmass):
        if s1_rockmass <= 0 or s1_induced <= 0:
            return 0
        else:
            return s1_rockmass/s1_induced


    # Función para plotear visualizaciones
    # No utilizada por la GUI, se deja para tener otra alternativa de gráfico
    def plot(self):
        columns = self.boundary_df.columns.values.tolist()
        columns.remove("x")
        columns.remove("y")

        print(columns)
        selection = input("Select variable to show: ")
        label_graphic = input("Write a label: ")
        min_value = input("Min value to legend: ")
        max_value = input("Max value to legend: ")

        graphic = plt.scatter(x = self.boundary_df["x"], 
                            y = self.boundary_df["y"], 
                            c = self.boundary_df[selection], 
                            cmap = "coolwarm", 
                            s = 1 #Tamaño de los puntos
                            ) 

        plt.colorbar(graphic, label = label_graphic, pad = 0.05)
        if min_value != "" and max_value != "":
            min_value = float(min_value)
            max_value = float(max_value)
            plt.clim(min_value, max_value)
        plt.grid(color='white', linestyle="--", linewidth=0.5)

        #Graficar geometrías circulares
        for geometry in self.geometries:
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
                    x_geom_upper[i] = self.transform_xt(x = x_geom[i]-geometry.x_centre, y = y_geom_u[i]-geometry.y_centre, theta = geometry.beta) + geometry.x_centre
                    y_geom_upper_contour[i] = self.transform_yt(x = x_geom[i]-geometry.x_centre, y = y_geom_u[i]-geometry.y_centre, theta = geometry.beta) + geometry.y_centre
                    x_geom_lower[i] = self.transform_xt(x = x_geom[i]-geometry.x_centre, y = y_geom_l[i]-geometry.y_centre, theta = geometry.beta) + geometry.x_centre
                    y_geom_lower_contour[i] = self.transform_yt(x = x_geom[i]-geometry.x_centre, y = y_geom_l[i]-geometry.y_centre, theta = geometry.beta) + geometry.y_centre


            plt.plot(x_geom_upper, y_geom_upper_contour, "k")
            plt.plot(x_geom_lower, y_geom_lower_contour, "k")
        plt.show()

    # Cálculo de esfuerzos utilizando fórmulas de Kirsch y guardado en dataframe
    def calculate_kirsch_stresses(self, P, Q, angle_P=0):
        self.P = P
        self.Q = Q
        self.angle_P = np.radians(angle_P) #Ángulo = 0° si es paralelo al eje X

        self.sx_insitu = self.rotate_sigma_x(sx = self.P, sy = self.Q, txy = 0, theta = -self.angle_P)
        self.sy_insitu = self.rotate_sigma_y(sx = self.P, sy = self.Q, txy = 0, theta = -self.angle_P)
        self.txy_insitu = self.rotate_tau_xy(sx = self.P, sy = self.Q, txy = 0, theta = -self.angle_P)

        #Cálculo de p y K para elipse
        self.p_ellipse = Q
        self.K_ellipse = P/Q

        for geometry in self.geometries:
            #En caso de ser circunferencia
            if geometry.__class__.__name__ == "Geometry_circle":
                #Listado de datos a guardar
                r = []
                theta = []
                sigma_r = []
                sigma_theta = []
                tau_r_theta = []

                for i in self.boundary_df.index:
                    #self.kirsch_r(x_point, y_point, x_centre_geometry, y_centre_geometry)
                    calc_r = self.kirsch_r(x_point = self.boundary_df["x"][i], 
                                        y_point = self.boundary_df["y"][i], 
                                        x_centre_geometry = geometry.x_centre, 
                                        y_centre_geometry = geometry.y_centre
                                        )
                    #self.kirsch_theta(x_point, y_point, x_centre_geometry, y_centre_geometry)
                    calc_theta = self.kirsch_theta(x_point = self.boundary_df["x"][i], 
                                                y_point = self.boundary_df["y"][i], 
                                                x_centre_geometry = geometry.x_centre, 
                                                y_centre_geometry = geometry.y_centre
                                                )
                    #self.kirsch_sigma_r(sx, sy, txy, r, a, theta)
                    calc_sigma_r = self.kirsch_sigma_r(sx = self.sx_insitu, 
                                                    sy = self.sy_insitu, 
                                                    txy = self.txy_insitu, 
                                                    r = calc_r, 
                                                    a = geometry.radius, 
                                                    theta = calc_theta
                                                    )
                    #self.kirsch_sigma_theta(sx, sy, txy, r, a, theta)
                    calc_sigma_theta = self.kirsch_sigma_theta(sx = self.sx_insitu, 
                                                            sy = self.sy_insitu, 
                                                            txy = self.txy_insitu, 
                                                            r = calc_r, 
                                                            a = geometry.radius, 
                                                            theta = calc_theta
                                                            )
                    #self.kirsch_tau_r_theta(sx, sy, txy, r, a, theta)
                    calc_tau_r_theta = self.kirsch_tau_r_theta(sx = self.sx_insitu, 
                                                            sy = self.sy_insitu, 
                                                            txy = self.txy_insitu, 
                                                            r = calc_r, 
                                                            a = geometry.radius, 
                                                            theta = calc_theta
                                                            )

                    #Guardado de datos por iteración en las listas temporales
                    r.append(calc_r)
                    theta.append(calc_theta)
                    sigma_r.append(calc_sigma_r)
                    sigma_theta.append(calc_sigma_theta)
                    tau_r_theta.append(calc_tau_r_theta)


                #Guardado de listas en la base de datos de la BoundaryBox    
                self.boundary_df[geometry.name+"_r"] = r
                self.boundary_df[geometry.name+"_theta"] = theta
                self.boundary_df[geometry.name+"_sigma_r"] = sigma_r
                self.boundary_df[geometry.name+"_sigma_theta"] = sigma_theta
                self.boundary_df[geometry.name+"_tau_r_theta"] = tau_r_theta

            elif geometry.__class__.__name__ == "Geometry_ellipse":
                #Listado de datos a guardar
                ###########Guardados auxiliares##############
                # Parámetros
                #e0 = []
                #b = []
                #d = []
                #u = []
                #e = []
                #phi = []
                theta = [] #Este si se debe guardar
                #C = []
                #J = []
                sigma_r = [] #Este si se debe guardar
                sigma_theta = [] #Este si se debe guardar
                tau_r_theta = [] #Este si se debe guardar

                for i in self.boundary_df.index:
                    #Transformación de coordenadas (x,y) -> (x',y')
                    x_t = self.transform_x(x = self.boundary_df["x"][i]-geometry.x_centre, y = self.boundary_df["y"][i]-geometry.y_centre, theta = self.angle_P) 
                    y_t = self.transform_y(x = self.boundary_df["x"][i]-geometry.x_centre, y = self.boundary_df["y"][i]-geometry.y_centre, theta = self.angle_P) 

                    #kirsch_ellipse_e0(self, W, H):
                    calc_e0 = self.kirsch_ellipse_e0(W = geometry.W, H = geometry.H)

                    #kirsch_ellipse_b(self, W, H, x1, z1):
                    calc_b = self.kirsch_ellipse_b(W = geometry.W, H = geometry.H, x1 = x_t, z1 = y_t)

                    #kirsch_ellipse_d(self, W, H, x1, z1):
                    calc_d = self.kirsch_ellipse_d(W = geometry.W, H = geometry.H, x1 = x_t, z1 = y_t)

                    #kirsch_ellipse_u(self, b, e0, d):
                    calc_u = self.kirsch_ellipse_u(b = calc_b, e0 = calc_e0, d = calc_d)

                    #kirsch_ellipse_e(self, e0, u):
                    calc_e = self.kirsch_ellipse_e(e0 = calc_e0, u = calc_u)

                    #kirsch_ellipse_phi(self, e, x1, z1):
                    calc_phi = self.kirsch_ellipse_phi(e = calc_e, x1 = x_t, z1 = y_t)

                    #kirsch_ellipse_theta(self, e, x1, z1):
                    calc_theta = self.kirsch_ellipse_theta(e = calc_e, x1 = x_t, z1 = y_t)

                    #kirsch_ellipse_C(self, e0, e):
                    calc_C = self.kirsch_ellipse_C(e0 = calc_e0, e = calc_e)

                    #kirsch_ellipse_J(self, e, phi):
                    calc_J = self.kirsch_ellipse_J(e = calc_e, phi = calc_phi)

                    # Cálculo de esfuerzos 
                    # "tf" se refiere a que calcula esfuerzos en la coordenada (x',y') transformada
                    #kirsch_ellipse_sigma_r(self, p, e0, e, C, J, K, phi, beta):
                    calc_sigma_r_tf = self.kirsch_ellipse_sigma_r(p = self.p_ellipse, e0 = calc_e0, e = calc_e, 
                                                                C = calc_C, J = calc_J, K = self.K_ellipse, 
                                                                phi = calc_phi, beta = geometry.beta
                                                                )
                    #kirsch_ellipse_sigma_theta(self, p, e0, e, J, K, phi, beta, sigma_r):
                    calc_sigma_theta_tf = self.kirsch_ellipse_sigma_theta(p = self.p_ellipse, e0 = calc_e0, e = calc_e,
                                                                        J = calc_J, K = self.K_ellipse, phi = calc_phi,
                                                                        beta = geometry.beta, sigma_r = calc_sigma_r_tf
                                                                        )
                    #kirsch_ellipse_tau_r_theta(self, p, e0, e, C, J, K, phi, beta):
                    calc_tau_r_theta_tf = self.kirsch_ellipse_tau_r_theta(p = self.p_ellipse, e0 = calc_e0, e = calc_e, 
                                                                        C = calc_C, J = calc_J, K = self.K_ellipse, 
                                                                        phi = calc_phi, beta = geometry.beta
                                                                        )
                    ###Rotación a sistema (x,y) (Rotación en base a inclinación de esfuerzo P -> angle_P)
                    calc_sigma_r = self.rotate_sigma_x(sx = calc_sigma_r_tf, 
                                                    sy = calc_sigma_theta_tf, 
                                                    txy = calc_tau_r_theta_tf, 
                                                    theta = -self.angle_P
                                                    )
                    calc_sigma_theta = self.rotate_sigma_y(sx = calc_sigma_r_tf, 
                                                        sy = calc_sigma_theta_tf, 
                                                        txy = calc_tau_r_theta_tf, 
                                                        theta = -self.angle_P
                                                        )
                    calc_tau_r_theta = self.rotate_tau_xy(sx = calc_sigma_r_tf, 
                                                        sy = calc_sigma_theta_tf, 
                                                        txy = calc_tau_r_theta_tf, 
                                                        theta = -self.angle_P
                                                        )
                    #Guardado de datos por iteración en las listas temporales
                    ######GUARDADOS AUXILIARES
                    #e0.append(calc_e0)
                    #b.append(calc_b)
                    #d.append(calc_d)
                    #u.append(calc_u)
                    #e.append(calc_e)
                    #phi.append(calc_phi)
                    #C.append(calc_C)
                    #J.append(calc_J)
                    ######GUARDADOS AUXILIARES

                    #Rotación para obtener Sx, Sy y Txy es beta + theta en el caso de elipse
                    theta.append(calc_theta + geometry.beta) 
                    sigma_r.append(calc_sigma_r)
                    sigma_theta.append(calc_sigma_theta)
                    tau_r_theta.append(calc_tau_r_theta)

                ######GUARDADOS AUXILIARES
                #self.boundary_df[geometry.name+"_e0"] = e0
                #self.boundary_df[geometry.name+"_b"] = b
                #self.boundary_df[geometry.name+"_d"] = d
                #self.boundary_df[geometry.name+"_u"] = u
                #self.boundary_df[geometry.name+"_e"] = e
                #self.boundary_df[geometry.name+"_phi"] = phi
                #self.boundary_df[geometry.name+"_C"] = C
                #self.boundary_df[geometry.name+"_J"] = J
                ######GUARDADOS AUXILIARES

                #Guardado de listas en la base de datos de la BoundaryBox 
                self.boundary_df[geometry.name+"_theta"] = theta
                self.boundary_df[geometry.name+"_sigma_r"] = sigma_r
                self.boundary_df[geometry.name+"_sigma_theta"] = sigma_theta
                self.boundary_df[geometry.name+"_tau_r_theta"] = tau_r_theta

        print("Kirsch's stresses calculated")

    #Cálculo de superposición de esfuerzos
    def calculate_stress_overlap(self):
        sx_t = []
        sy_t = []
        txy_t = []
        for i in self.boundary_df.index:
            calc_sx = self.sx_insitu
            calc_sy = self.sy_insitu
            calc_txy = self.txy_insitu
            
            #Rotación de esfuerzos por geometría
            for geometry in self.geometries:
                rotated_sx_geometry = self.rotate_sigma_x(sx = self.boundary_df[geometry.name+"_sigma_r"][i],
                                                        sy = self.boundary_df[geometry.name+"_sigma_theta"][i],
                                                        txy = self.boundary_df[geometry.name+"_tau_r_theta"][i],
                                                        theta = -self.boundary_df[geometry.name+"_theta"][i]
                                                        )
                rotated_sy_geometry = self.rotate_sigma_y(sx = self.boundary_df[geometry.name+"_sigma_r"][i],
                                                        sy = self.boundary_df[geometry.name+"_sigma_theta"][i],
                                                        txy = self.boundary_df[geometry.name+"_tau_r_theta"][i],
                                                        theta = -self.boundary_df[geometry.name+"_theta"][i]
                                                        )
                rotated_txy_geometry = self.rotate_tau_xy(sx = self.boundary_df[geometry.name+"_sigma_r"][i],
                                                        sy = self.boundary_df[geometry.name+"_sigma_theta"][i],
                                                        txy = self.boundary_df[geometry.name+"_tau_r_theta"][i],
                                                        theta = -self.boundary_df[geometry.name+"_theta"][i]
                                                        )

                #Suma de contribuciones por geometría en cada componente
                calc_sx += rotated_sx_geometry - self.sx_insitu
                calc_sy += rotated_sy_geometry - self.sy_insitu
                calc_txy += rotated_txy_geometry - self.txy_insitu

            #Guardado de datos en listas temporales
            sx_t.append(calc_sx)
            sy_t.append(calc_sy)
            txy_t.append(calc_txy)

        #Guardado de cálculos en dataframe
        self.boundary_df["Sx"] = sx_t
        self.boundary_df["Sy"] = sy_t
        self.boundary_df["Txy"] = txy_t

        print("Sx, Sy and Txy calculated")

    # Cálculo de esfuerzos principales
    def calculate_principal_stresses(self):
        stress_P = []
        stress_Q = []
        stress_max_shear = []

        for i in self.boundary_df.index:
            calc_p = self.principal_stress_P(sx = self.boundary_df["Sx"][i], 
                                            sy = self.boundary_df["Sy"][i], 
                                            txy = self.boundary_df["Txy"][i]
                                            )
            calc_q = self.principal_stress_Q(sx = self.boundary_df["Sx"][i], 
                                            sy = self.boundary_df["Sy"][i], 
                                            txy = self.boundary_df["Txy"][i]
                                            )
            calc_tmax = self.principal_stress_max_shear(sx = self.boundary_df["Sx"][i], 
                                                        sy = self.boundary_df["Sy"][i], 
                                                        txy = self.boundary_df["Txy"][i]
                                                        )

            stress_P.append(calc_p)
            stress_Q.append(calc_q)
            stress_max_shear.append(calc_tmax)

        self.boundary_df["P"] = stress_P
        self.boundary_df["Q"] = stress_Q
        self.boundary_df["Tmax"] = stress_max_shear

        print("Principal stresses calculated")

    # Cálculo de factores de seguridad
    def calculate_strength_factor(self, sci, mi, gsi):
        P_rockmass = []
        strength_factor = []

        for i in self.boundary_df.index:
            #hoek_brown_S1(self, s3, sci_rock, mi_rock, gsi_rockmass)
            calc_p_rockmass = self.hoek_brown_S1(s3 = self.boundary_df["Q"][i],
                                                sci_rock = sci,
                                                mi_rock = mi,
                                                gsi_rockmass = gsi
                                                )
            #strength_factor(self, s1_induced, s1_rockmass)
            calc_strength_factor = self.strength_factor(s1_induced = self.boundary_df["P"][i], 
                                                        s1_rockmass = calc_p_rockmass
                                                        )
            P_rockmass.append(calc_p_rockmass)
            strength_factor.append(calc_strength_factor)

        self.boundary_df["P_rockmass"] = P_rockmass
        self.boundary_df["Strength_factor"] = strength_factor

        print("Strength factor calculated")

    # Método para exportar dataframe a .csv
    # La GUI no utiliza este método, pero está a modo de alternativa
    def export_to_csv(self, name_file):
        #Añadir en variables de entrada a "parameters" = []

        #Ver como incorporar la ruta del archivo
        self.boundary_df.to_csv(name_file + ".csv", sep = ",")
        print("Export successfully")

# Ejemplo de uso de backend
############# Prueba de Cálculos
#A = BoundaryBox(70.5, 70.5, 0.5, 0.5)
#B = Geometry_circle("B", 30.5, 30.5, 5.5)
#C = Geometry_ellipse("C", -10.3, -10.3, 15.2, 7.7, 30)
#A.mesh_points([B, C])
#A.calculate_kirsch_stresses(30, 10, 20)
#A.calculate_stress_overlap()
#A.calculate_principal_stresses()
#A.plot()
