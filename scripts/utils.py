import os
import csv
import numpy as np
import math
from scipy.interpolate import interp1d


class Utils:
    """
    Utility class for processing aerodynamic data, performing numerical analysis, 
    and computing stability and flight quality metrics.
    """

    def __init__(self, file_number):
        """
        Initialize the Utils instance with a file identifier.

        Args:
            file_number (int): Identifier for the file being processed.
        """
        self.file_number = file_number

    def save_lists_from_input(self, list_of_lists, input_strings):
        """
        Save sublists that contain all the provided input strings until an empty list is encountered.

        Args:
            list_of_lists (list of lists): List containing sublists to search through.
            input_strings (iterable): Strings to search for in each sublist.

        Returns:
            list: A list of sublists that meet the criteria.
        """
        saved_lists = []  # Here we'll store the lists that meet the criteria
        found_input = False

        for sublist in list_of_lists:
            if all(input_string in sublist for input_string in input_strings):
                found_input = True

            if found_input:
                saved_lists.append(sublist)

            # Stop if an empty list is found after the criteria is met
            if found_input and not sublist:
                break
        return saved_lists

    def extract_columns(self, data):
        """
        Extract numeric columns from the provided data. Determines whether to extract 2 or 6 columns.

        Args:
            data (list of lists): Data rows to extract columns from.

        Returns:
            list: A list containing the extracted numeric columns.
        """
        if len(data[1]) != 6:
            num_columns = 2 
        else:
            num_columns = 6 
        columns = [[] for _ in range(num_columns)]
        for row in data:
            if len(row) == num_columns:
                try:
                    for i in range(num_columns):
                        if row[i] != "":
                            columns[i].append(float(row[i]))
                        else:
                            columns[i].append(0.0)
                except ValueError:
                    pass

        return columns

    def search_in_out(self, lines, data_dict, alpha_cruise):
        """
        Search through lines from an output file to extract aerodynamic parameters,
        interpolate data, and update the provided data dictionary.

        Args:
            lines (list): Lines from the output file.
            data_dict (dict): Dictionary containing previously extracted data.
            alpha_cruise (float): Angle of attack for cruise conditions used in interpolation.

        Returns:
            dict: Updated data dictionary with extracted and interpolated parameters.
        """
        count_zero_l_a = 0
        count_config = 0
        alpha_len = len(data_dict["Alpha_list"])

        for idx, line in enumerate(lines):
            if "ZERO LIFT ANGLE OF ATTACK" in line:
                list_aux = line.split(" ")
                if count_zero_l_a == 0:
                    data_dict["Zero_lift_angle_w"] = [float(i.replace("\n", "")) for i in list_aux 
                                                      if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]

                elif count_zero_l_a == 1:
                    data_dict["Zero_lift_angle_htp"] = [float(i.replace("\n", "")) for i in list_aux 
                                                        if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]

                elif count_zero_l_a == 2:
                    data_dict["Zero_lift_angle_vtp"] = [float(i.replace("\n", "")) for i in list_aux 
                                                        if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]

            elif "ZERO LIFT PITCHING MOMENT COEFFICIENT" in line:
                list_aux = line.split(" ")
                if count_zero_l_a == 0:
                    data_dict["Zero_lift_pmc_w"] = [float(i.replace("\n", "")) for i in list_aux 
                                                    if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]
                    count_zero_l_a += 1

                elif count_zero_l_a == 1:
                    data_dict["Zero_lift_pmc_htp"] = [float(i.replace("\n", "")) for i in list_aux 
                                                      if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]
                    count_zero_l_a += 1

                elif count_zero_l_a == 2:
                    data_dict["Zero_lift_pmc_vtp"] = [float(i.replace("\n", "")) for i in list_aux 
                                                      if i.replace("-", "").replace(".", "").replace("E", "").replace("\n", "").isdigit()][0]
                    count_zero_l_a += 1

            elif ("CHARACTERISTICS AT ANGLE OF ATTACK AND IN SIDESLIP" in line and 
                  "WING ALONE CONFIGURATION" in lines[idx+1] and count_config == 0):
                list_aux = lines[idx:]
                matrix = self.process_text_matrix(list_aux, alpha_len)
                awb = self.interpolate_fun(data_dict, "Alpha_list", matrix[7], alpha_cruise)

                data_dict["awb"] = float(awb) * 3.1415 / 180
                count_config += 1

                list_aux = lines[idx+8].split(" ")
                count_col = 0
                for i in list_aux:
                    try:
                        data_dict["MAC"] = float(i)
                        count_col += 1
                        if count_col == 9:
                            break
                    except:
                        pass
                count_col = 0
                for i in list_aux:
                    try:
                        data_dict["Sref"] = float(i)
                        count_col += 1
                        if count_col == 8:
                            break
                    except:
                        pass

            elif ("CHARACTERISTICS AT ANGLE OF ATTACK AND IN SIDESLIP" in line and 
                  "HORIZONTAL TAIL CONFIGURATION" in lines[idx+1] and count_config == 1):
                list_aux = lines[idx:]
                matrix = self.process_text_matrix(list_aux, alpha_len)
                atStS = self.interpolate_fun(data_dict, "Alpha_list", matrix[7], alpha_cruise)

                data_dict["at*St/S"] = float(atStS) * 3.1415 / 180
                count_config += 1

        return data_dict

    def process_text_matrix(self, text, alpha_len):
        """
        Process a block of text lines to form a transposed numeric matrix.

        Args:
            text (list): Lines of text to be processed.
            alpha_len (int): Number of rows to extract based on alpha length.

        Returns:
            list: Transposed matrix with numeric values.
        """
        # Extract the relevant lines from the text block
        lines = text[12:12+alpha_len]

        # Split each line into elements
        matrix = [line.split() for line in lines]

        # Convert elements to floats, or 0.0 if conversion fails
        matrix_floats = []
        for row in matrix:
            row_floats = []
            for element in row:
                try:
                    row_floats.append(float(element))
                except ValueError:
                    row_floats.append(0.0)
            matrix_floats.append(row_floats)

        # Transpose the matrix
        matrix_floats_transposed = list(map(list, zip(*matrix_floats)))
        return matrix_floats_transposed

    def read_from_csv(self, csv_path, flight_conditions, to_trim, geometry_dict, alpha_cruise, mass, prev_data_dict):
        """
        Read aerodynamic data from a CSV file and a corresponding output file, 
        process the data, perform interpolations, and compute trim values.

        Args:
            csv_path (str): Path to the CSV file.
            flight_conditions (tuple): Flight conditions (e.g., altitude, density, velocity).
            to_trim (any): Indicator for whether trim calculation is needed.
            geometry_dict (dict): Dictionary containing aircraft geometry.
            alpha_cruise (float): Angle of attack during cruise for interpolation.
            mass (float): Aircraft mass.
            prev_data_dict (dict): Previous data dictionary to fall back on if necessary.

        Returns:
            dict: Data dictionary containing processed aerodynamic parameters.
        """
        column_names = [
            ("Alpha", "Cd"), ("Alpha", "CL"), ("Alpha", "Cm"), ("Alpha", "Cyb"), ("Alpha", "Cnb"), ("Alpha", "Clb"),
            ("Alpha", "Q_Qinf"), ("Alpha", "Epslon"), ("Alpha", "d_Epslon"), ("Alpha", "CLq"), ("Alpha", "Cmq"),
            ("Alpha", "CLad_Basic"), ("Alpha", "Cmad_Basic"), ("Alpha", "Clp"), ("Alpha", "Cyp"), ("Alpha", "Cnp"),
            ("Alpha", "Cnr"), ("Alpha", "Clr"), ("DeltaFlap", "CL"), ("DeltaFlap", "Cm"),
            ("DeltaElev", "CL_Elev"), ("DeltaElev", "Cm_Elev"), ("DeltaElev", "CLa_Elev"),
            ("DeltaElev", "Cha_Elev"), ("DeltaElev", "Chd_Elev")
        ]

        data_dict = {}
        # Read CSV file and extract data lists
        with open(csv_path, newline='') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',', skipinitialspace=True)
            data_lists = [row for row in csvreader]
            for column_name in column_names:
                column_data = self.extract_columns(self.save_lists_from_input(data_lists, column_name))
                if len(column_data) == 2:
                    data_dict[column_name[0] + "_list"] = column_data[0]
                    data_dict[column_name[0] + "_" + column_name[1] + "_list"] = column_data[1]

            data = self.extract_columns(self.save_lists_from_input(data_lists, ("Alpha", "DeltaFlap")))
            data_dict["Alpha_DeltaFlap_dict"] = {
                data[1][0]: data[1][1:], 
                data[2][0]: data[2][1:], 
                data[3][0]: data[3][1:], 
                data[4][0]: data[4][1:], 
                data[5][0]: data[5][1:]
            }
            data = self.extract_columns(self.save_lists_from_input(data_lists, ["Alpha", "DeltaElev"]))
            data_dict["Alpha_DeltaElev_dict"] = {
                data[1][0]: data[1][1:], 
                data[2][0]: data[2][1:], 
                data[3][0]: data[3][1:], 
                data[4][0]: data[4][1:], 
                data[5][0]: data[5][1:]
            }
        with open(csv_path.replace(".csv", ".out"), "r") as outfile:
            lines = outfile.readlines()

        data_dict = self.search_in_out(lines, data_dict, alpha_cruise)

        # Fallback to previous data if incomplete
        if len(data_dict) != 40:
            data_dict = prev_data_dict

        # Calculate trim if required
        if to_trim is None:
            data_dict["trim"] = [0.0, 0.0]

            zero_lift_pmc_w = data_dict["Zero_lift_pmc_w"]
            zero_lift_pmc_htp = data_dict["Zero_lift_pmc_htp"]
            mass = mass
            xcg = geometry_dict["XCG_GEN"]
            rho = flight_conditions[2]
            velocity = flight_conditions[3]
            Sref = data_dict["Sref"]
            atStS = data_dict["at*St/S"]
            awb = data_dict["awb"]
            zero_lift_angle_w = data_dict["Zero_lift_angle_w"]
            zero_lift_angle_htp = data_dict["Zero_lift_angle_htp"]
            cw = 2 * mass * 9.81 / (rho * (velocity ** 2) * Sref)
            mac = data_dict["MAC"]
            xwb = geometry_dict["XW_GEN"]
            xt = geometry_dict["XH_GEN"]
            xwbbar = (xcg - xwb) / mac
            xtbar = (xt - xcg) / mac

            Qinf = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Q_Qinf_list", alpha_cruise)
            if Qinf == float("nan") or math.isnan(Qinf):
                Qinf = 1.0

            epsilon = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Epslon_list", alpha_cruise)
            if epsilon == float("nan") or math.isnan(epsilon):
                epsilon = 0.0

            # Calculate horizontal tail incidence (i_t)
            it = (zero_lift_pmc_w + zero_lift_pmc_htp + cw * xwbbar) / (atStS * Qinf * ((xt - xwb) / mac)) - alpha_cruise + epsilon + zero_lift_angle_htp
            data_dict["trim"][1] = self.round_to_significant_figures(it, 6)
            
            # Calculate wing-body incidence (i_wb)
            iwb = (cw - atStS * Qinf * (alpha_cruise - epsilon + it - zero_lift_angle_htp)) / awb - alpha_cruise + zero_lift_angle_w
            data_dict["trim"][0] = self.round_to_significant_figures(iwb, 6)
            
        return data_dict

    def str_to_val(self, data):
        """
        Convert all numerical values in the dictionary from strings or integers to floats.

        Args:
            data (dict): Dictionary with values that may be strings, integers, or lists.

        Returns:
            dict: Updated dictionary with all values converted to floats.
        """
        for key, value in data.items():
            if isinstance(value, (int, float)):
                data[key] = float(value)
            elif isinstance(value, list):
                data[key] = [float(num) for num in value]
            elif isinstance(value, str):
                data[key] = float(value)
        return data

    def interpolate_fun(self, data_dict, x, y, alpha_cruise):
        """
        Perform linear interpolation for the given x and y data based on alpha_cruise.

        Args:
            data_dict (dict): Dictionary containing the x data if x or y is a string.
            x (list or str): Independent variable or key for the independent variable.
            y (list or str): Dependent variable or key for the dependent variable.
            alpha_cruise (float): The interpolation point (angle of attack in cruise).

        Returns:
            float: Interpolated value at alpha_cruise.
        """
        if isinstance(x, list) and isinstance(y, list):
            return np.interp(alpha_cruise, x, y)
        elif isinstance(x, str) and isinstance(y, list):
            return np.interp(alpha_cruise, data_dict[x], y)
        elif isinstance(x, list) and isinstance(y, str):
            return np.interp(alpha_cruise, x, data_dict[y])
        else:
            return np.interp(alpha_cruise, data_dict[x], data_dict[y])

    def calculate_derivative(self, data_dict, x, y, point):
        """
        Calculate the numerical derivative of y with respect to x at a given point using centered difference.

        Args:
            data_dict (dict): Dictionary containing x and y data lists.
            x (str): Key for the independent variable in data_dict.
            y (str): Key for the dependent variable in data_dict.
            point (float): The point at which the derivative is computed.

        Returns:
            float: Approximated derivative value at the given point.
        """
        x_data, y_data = data_dict[x], data_dict[y]
        # Interpolate the data using linear interpolation
        interpolation = interp1d(x_data, y_data, kind='linear')

        # Calculate values above and below the point
        delta_x = 0.01
        x_previous = point - delta_x
        x_next = point + delta_x

        y_previous = interpolation(x_previous)
        y_next = interpolation(x_next)

        try:
            derivative = (y_next - y_previous) / (x_next - x_previous)
        except Exception as e:
            print("ERROR DERIVATIVE: ", e)
            derivative = 0.0

        return derivative

    def round_to_significant_figures(self, number, significant_figures):
        """
        Round a number to the specified number of significant figures.

        Args:
            number (float): The number to be rounded.
            significant_figures (int): The number of significant figures desired.

        Returns:
            float: The rounded number.
        """
        try:
            if number == 0:
                return 0
            else:
                return round(number, significant_figures - int(math.floor(math.log10(abs(number)))) - 1)
        except:
            return 0.0

    def get_stability_derivatives(self, data_dict, geometry_dict, alpha_cruise, mach, x):
        """
        Compute stability derivatives using interpolated aerodynamic data and geometry parameters.

        Args:
            data_dict (dict): Dictionary containing aerodynamic data.
            geometry_dict (dict): Dictionary containing aircraft geometry parameters.
            alpha_cruise (float): Angle of attack during cruise for interpolation.
            mach (float): Mach number.
            x (list): List of additional geometric parameters.

        Returns:
            dict: Dictionary of computed stability derivatives, rounded to six significant figures.
        """
        sspn_w = x[4]
        x_cg = x[0]
        fuselage_len = x[1]
        derivatives_dict = {}
        derivatives_dict["Cxu"] = -2 * self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cd_list", alpha_cruise)
        derivatives_dict["CLs"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_CL_list", alpha_cruise)
        derivatives_dict["Czs"] = -derivatives_dict["CLs"]
        derivatives_dict["Cds"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cd_list", alpha_cruise)
        derivatives_dict["Czu"] = (-mach**2 / (1 - mach**2)) * derivatives_dict["CLs"]
        derivatives_dict["Cmu"] = 0.0
        derivatives_dict["CLa"] = self.calculate_derivative(data_dict, "Alpha_list", "Alpha_CL_list", alpha_cruise) * 180 / np.pi
        derivatives_dict["Cda"] = self.calculate_derivative(data_dict, "Alpha_list", "Alpha_Cd_list", alpha_cruise) * 180 / np.pi
        derivatives_dict["Cxa"] = -derivatives_dict["Cda"] + derivatives_dict["CLs"]
        derivatives_dict["Cza"] = -derivatives_dict["CLa"] - derivatives_dict["Cds"]
        derivatives_dict["Cma"] = self.calculate_derivative(data_dict, "Alpha_list", "Alpha_Cm_list", alpha_cruise) * 180 / np.pi
        derivatives_dict["CLq"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_CLq_list", alpha_cruise)
        derivatives_dict["Czq"] = -derivatives_dict["CLq"]
        derivatives_dict["Cmq"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cmq_list", alpha_cruise)
        derivatives_dict["CLad"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_CLad_Basic_list", alpha_cruise)
        derivatives_dict["Cmad"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cmad_Basic_list", alpha_cruise)
        derivatives_dict["Czad"] = -derivatives_dict["CLad"]
        derivatives_dict["Cyb"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cyb_list", alpha_cruise)
        derivatives_dict["Clb"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Clb_list", alpha_cruise)
        derivatives_dict["Cnb"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cnb_list", alpha_cruise)
        derivatives_dict["Cyp"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cyp_list", alpha_cruise)
        derivatives_dict["Clp"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Clp_list", alpha_cruise)
        derivatives_dict["Cnp"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cnp_list", alpha_cruise)
        derivatives_dict["Cnr"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Cnr_list", alpha_cruise)
        derivatives_dict["Cyr"] = -0.5 * derivatives_dict["Cnr"] / ((0.97 * fuselage_len - x_cg) / (2 * sspn_w)) 
        derivatives_dict["Clr"] = self.interpolate_fun(data_dict, "Alpha_list", "Alpha_Clr_list", alpha_cruise)

        for key, value in derivatives_dict.items():
            derivatives_dict[key] = self.round_to_significant_figures(value, 6)

        return derivatives_dict

    def get_flight_quality(self, data_dict, derivatives_dict, geometry_dict, flight_conditions, mass, semi_span_w):
        """
        Compute flight quality metrics for both longitudinal and lateral-directional dynamics.

        Args:
            data_dict (dict): Dictionary containing aerodynamic parameters.
            derivatives_dict (dict): Dictionary containing stability derivatives.
            geometry_dict (dict): Dictionary containing aircraft geometry.
            flight_conditions (tuple): Tuple containing flight conditions (e.g., altitude, density, velocity).
            mass (float): Aircraft mass.
            semi_span_w (float): Wing semi-span.

        Returns:
            dict: Dictionary containing flight quality metrics for longitudinal (CP, FUG) and lateral-directional (SP, CB, BH) modes.
        """
        flight_quality_dict = {}
        # Longitudinal calculations
        m = mass
        fuselage_len = geometry_dict["X(1)_BODY"][-1]
        c = data_dict["MAC"]
        g = 9.81
        rho = flight_conditions[2]
        S = data_dict["Sref"]
        mu_long = m / (0.5 * rho * S * c)
        Iyy = 6205.5  # Moment of inertia about the y-axis
        V = flight_conditions[3]
        t_ref = c / (2 * V)

        CZa_dot = derivatives_dict["Czad"]
        Cma_dot = derivatives_dict["Cmad"]
        CXu = derivatives_dict["Cxu"]
        CXa = derivatives_dict["Cxa"]
        CZs = derivatives_dict["Czs"]
        CZu = derivatives_dict["Czu"]
        CZa = derivatives_dict["Cza"]
        CZq = derivatives_dict["Czq"]
        Cma = derivatives_dict["Cma"]
        Cmq = derivatives_dict["Cmq"]
        CLa = derivatives_dict["CLa"]

        M_long = np.array([[2 * mu_long, 0, 0, 0],
                             [0, 2 * mu_long - CZa_dot, 0, 0],
                             [0, -Cma_dot, Iyy, 0],
                             [0, 0, 0, 1]])

        A_long = np.array([[CXu, CXa, 0, CZs],
                             [2 * CZs + CZu, CZa, 2 * mu_long + CZq, 0],
                             [0, Cma, Cmq, 0],
                             [0, 0, 1, 0]])

        A_long = np.linalg.solve(M_long, A_long)

        autoval, autovect = np.linalg.eig(A_long)

        autoval = autoval / t_ref
        # Short period mode
        n1 = np.real(autoval[0])
        w1 = np.abs(np.imag(autoval[0]))
        w_n1 = np.sqrt(n1 ** 2 + w1 ** 2)
        amort1 = -n1 / w_n1
        try:
            T1 = 2 * np.pi / w1
        except:
            T1 = -1
        t_mitad1 = np.log(1 / 2) / n1

        nalpha = 0.5 * rho * V ** 2 * S * CLa / (m * g)
        if nalpha < 0:
            nalpha = 0
        wnspmax = np.sqrt(3.6 * nalpha)
        wnspmin = np.sqrt(0.085 * nalpha)

        # Phugoid mode
        n2 = np.real(autoval[3])
        w2 = np.abs(np.imag(autoval[3]))
        w_n2 = np.sqrt(n2 ** 2 + w2 ** 2)
        amort2 = -n2 / w_n2
        try:
            T2 = 2 * np.pi / w2
        except:
            T2 = -1
        t_mitad2 = np.log(1 / 2) / n2

        flight_quality_dict["CP"] = {
            "w_n": self.round_to_significant_figures(w_n1, 4),
            "amort": self.round_to_significant_figures(amort1, 4),
            "t_mitad": self.round_to_significant_figures(t_mitad1, 4),
            "nalpha": self.round_to_significant_figures(nalpha, 4),
            "wnspmin": self.round_to_significant_figures(wnspmin, 4),
            "wnspmax": self.round_to_significant_figures(wnspmax, 4),
            "T": self.round_to_significant_figures(T1, 4)
        }

        flight_quality_dict["FUG"] = {
            "amort": self.round_to_significant_figures(amort2, 4),
            "w_n": self.round_to_significant_figures(w_n2, 4),
            "t_mitad": self.round_to_significant_figures(t_mitad2, 4),
            "T": self.round_to_significant_figures(T2, 4)
        }

        # Lateral-Directional calculations
        b = semi_span_w * 2
        mu_latdir = 2 * m / (rho * S * b)
        Ix = 4.7970
        Ixz = -0.33281
        Iy = 12.18427
        Iz = 15.06720

        Czs = derivatives_dict["Czs"]
        Cyb = derivatives_dict["Cyb"]
        Cyp = derivatives_dict["Cyp"]
        Cyr = derivatives_dict["Cyr"]
        Clb = derivatives_dict["Clb"]
        Clp = derivatives_dict["Clp"]
        Clr = derivatives_dict["Clr"]
        Cnb = derivatives_dict["Cnb"]
        Cnp = derivatives_dict["Cnp"]
        Cnr = derivatives_dict["Cnr"]

        matriz_masa_latdir = np.array([[2 * mu_latdir, 0, 0, 0],
                                         [0, Ix, -Ixz, 0],
                                         [0, -Ixz, Iz, 0],
                                         [0, 0, 0, 1]])

        A_latdir = np.array([[Cyb, Cyp, Cyr - 2 * mu_latdir, -Czs],
                              [Clb, Clp, Clr, 0],
                              [Cnb, Cnp, Cnr, 0],
                              [0, 1, 0, 0]])

        A_latdir = np.linalg.solve(matriz_masa_latdir, A_latdir)

        autoval, autovectores = np.linalg.eig(A_latdir)
        autoval = autoval * (2 * V / b)

        reales = []
        complejos = []

        for i in range(4):
            if np.isreal(autoval[i]):
                reales.append(autoval[i])
            else:
                complejos.append(autoval[i])

        if len(reales) > 0:
            lambdaME = min(np.abs(reales))
            tMitadME = -np.log(1 / 2) / lambdaME
            tDobleME = -np.log(2) / lambdaME
            lambdaCB = max(np.abs(reales))
            tMitadCB = -np.log(1 / 2) / lambdaCB
        else:
            lambdaME = 0
            tMitadME = np.inf
            tDobleME = np.inf
            lambdaCB = 0
            tMitadCB = np.inf

        if len(complejos) > 0:
            nBH = np.real(complejos[0])
            wBH = np.abs(np.imag(complejos[0]))
            wnBH = np.sqrt(nBH ** 2 + wBH ** 2)
            amortBH = -nBH / wnBH
            T3 = 2 * np.pi / wBH
            tmitadBH = np.log(1 / 2) / nBH
        else:
            nBH = 0
            wBH = 0
            wnBH = 0
            amortBH = 0
            T3 = np.inf
            tmitadBH = np.inf

        flight_quality_dict["SP"] = {
            "t_doble": self.round_to_significant_figures(tDobleME, 4),
            "t_mitad": self.round_to_significant_figures(tMitadME, 4)
        }
        flight_quality_dict["CB"] = {
            "lambda": self.round_to_significant_figures(lambdaCB, 4),
            "t_mitad": self.round_to_significant_figures(tMitadCB, 4)
        }
        flight_quality_dict["BH"] = {
            "w_n": self.round_to_significant_figures(wnBH, 4),
            "amort": self.round_to_significant_figures(amortBH, 4),
            "t_mitad": self.round_to_significant_figures(tmitadBH, 4),
            "T": self.round_to_significant_figures(T3, 4)
        }

        return flight_quality_dict
