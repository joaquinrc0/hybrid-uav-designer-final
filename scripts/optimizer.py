"""
This module performs aerodynamic optimization using DATCOM.
It includes a loop to execute DATCOM, read simulation data, and an optimization process
that adjusts aircraft geometry to meet flight quality constraints.
"""

import os
import subprocess
import time

from scipy.optimize import minimize
import numpy as np

from scripts.utils import Utils
from scripts.dcm_creator_long import DCMCreator


class LoopThread:
    """
    Class to handle running the DATCOM simulation in a loop.
    
    This class constructs the necessary file paths based on a file number and
    provides methods to run DATCOM with given inputs and retrieve simulation data.
    
    Attributes:
        file_number (int): The identifier for file selection.
        methods (Utils): An instance of the Utils class for common utility functions.
        datcom_path (str): Path to the DATCOM executable batch file.
        dcm_path (str): Path to the DCM file used as input for DATCOM.
        csv_path (str): Path to the CSV file generated by DATCOM.
        working_dir (str): Working directory where DATCOM is executed.
        geometry_dict (dict): Dictionary to store geometry parameters extracted from DCM.
    """

    def __init__(self):
        """
        Initialize the LoopThread instance by setting file paths and creating a Utils instance.
        """
        # Extract the file number from the current script's name (hardcoded here as 2)
        self.file_number = 2
        self.methods = Utils(self.file_number)
        self.datcom_path = os.path.join(r"C:\Users\joroc\Desktop\1-MIA\TFM\Scripts\DATCOM", str(self.file_number), "bin\datcom.bat")
        self.dcm_path = os.path.join(r"C:\Users\joroc\Desktop\1-MIA\TFM\Scripts\DATCOM", str(self.file_number), "Example" + str(self.file_number), "eje" + str(self.file_number) + ".dcm")
        self.csv_path = os.path.join(r"C:\Users\joroc\Desktop\1-MIA\TFM\Scripts\DATCOM", str(self.file_number), "Example" + str(self.file_number), "eje" + str(self.file_number) + ".csv")
        self.working_dir = os.path.join(r"C:\Users\joroc\Desktop\1-MIA\TFM\Scripts\DATCOM", str(self.file_number), "Example" + str(self.file_number))
        self.geometry_dict = {}

    def run_datcom(self, inputs, flight_conditions, alpha_cruise, mass, prev_data_dict, to_trim):
        """
        Run DATCOM with the specified inputs and retrieve simulation data.
        
        This method creates a DCM file based on the input parameters, executes DATCOM,
        and then reads the resulting CSV file to extract simulation data.
        
        Args:
            inputs (list): List of input parameters for DCM creation.
            flight_conditions (tuple): Flight condition parameters (e.g., Mach, altitude, etc.).
            alpha_cruise (float): Cruise angle of attack.
            mass (float): Aircraft mass.
            prev_data_dict (dict): Dictionary containing previous simulation data.
            to_trim: Parameter to control trimming; can be None or provided data.
        
        Returns:
            The data read from the CSV file after DATCOM execution.
        """
        create_dcm = DCMCreator(
            self.dcm_path,
            inputs[0], inputs[1], inputs[2], inputs[3], inputs[4],
            inputs[5], inputs[6], inputs[7], inputs[8], inputs[9],
            to_trim
        )
        create_dcm.write_dcm_file()
        self.geometry_dict = self.methods.str_to_val(create_dcm.geometry_dict)

        # Change directory to working directory for DATCOM execution
        os.chdir(self.working_dir)
        # Command to execute DATCOM and pass the DCM file as an argument
        command = [self.datcom_path, self.dcm_path]
        
        # Execute the command and capture stdout and stderr
        process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, universal_newlines=True
        )
        
        # Send a newline to continue DATCOM execution
        process.stdin.write("\n")
        process.stdin.flush()
        
        stdout, stderr = process.communicate()

        # Read data from the generated CSV file and return it
        return self.methods.read_from_csv(self.csv_path, flight_conditions, to_trim, self.geometry_dict, alpha_cruise, mass, prev_data_dict)

    def loop(self, inputs, flight_conditions, alpha_cruise, mass, prev_data_dict, to_trim):
        """
        Run the DATCOM simulation loop and return simulation data along with geometry.
        
        This method calls run_datcom and then returns the obtained data and the updated geometry.
        
        Args:
            inputs (list): List of input parameters for DATCOM.
            flight_conditions (tuple): Flight condition parameters.
            alpha_cruise (float): Cruise angle of attack.
            mass (float): Aircraft mass.
            prev_data_dict (dict): Previous simulation data.
            to_trim: Trimming parameter, can be None or specific data.
        
        Returns:
            tuple: A tuple containing the simulation data and geometry dictionary.
        """
        data = self.run_datcom(inputs, flight_conditions, alpha_cruise, mass, prev_data_dict, to_trim)
        return data, self.geometry_dict


class AerodynamicOptimization:
    """
    Class to perform aerodynamic optimization using SLSQP.
    
    This class sets up and runs the optimization process, adjusting aircraft geometry
    to optimize flight qualities based on simulation data from DATCOM.
    
    Attributes:
        thread_stop_event (threading.Event): Event to signal when to stop the optimization thread.
        x0 (list): Initial guess for optimization variables.
        initial_bounds (list): Original bounds for the optimization variables.
        bounds (list): Current bounds for the optimization variables.
        flight_conditions (tuple): Flight condition parameters.
        alpha_cruise (float): Cruise angle of attack.
        plot_class: Object responsible for plotting aircraft geometry during optimization.
        default_limits (dict): Dictionary of default limits for various flight parameters.
        mass (float): Aircraft mass.
        eps (float): Tolerance value for optimization.
        max_iter (int): Maximum number of iterations for the optimization.
        thread (LoopThread): An instance of LoopThread to run DATCOM.
        methods (Utils): Utility methods from the Utils class.
        prev_data_dict (dict): Stores previous simulation data.
    """

    def __init__(self, thread_stop_event, x0, bounds, plot_class, flight_conditions, alpha_cruise, main_var_entries, default_limits):
        """
        Initialize the AerodynamicOptimization instance with provided parameters.
        
        Args:
            thread_stop_event (threading.Event): Event to signal thread termination.
            x0 (list): Initial guess for the design variables.
            bounds (list): List of tuples specifying lower and upper bounds for each variable.
            plot_class: Plotting class instance to update the GUI.
            flight_conditions (tuple): Flight condition parameters (e.g., Mach, altitude, etc.).
            alpha_cruise (float): Cruise angle of attack.
            main_var_entries (list): List of Tkinter entry widgets containing flight parameters.
            default_limits (dict): Dictionary of default variable limits.
        """
        self.thread_stop_event = thread_stop_event
        self.x0 = x0
        self.initial_bounds = bounds
        self.bounds = bounds
        self.flight_conditions = flight_conditions
        self.alpha_cruise = alpha_cruise
        self.plot_class = plot_class
        self.default_limits = default_limits
        self.mass = float(main_var_entries[5].get())
        self.eps = float(main_var_entries[4].get())
        self.max_iter = int(main_var_entries[3].get())

        # Instantiate LoopThread and related utilities
        self.thread = LoopThread()
        self.methods = self.thread.methods
        file_number = self.thread.file_number
        data_saving_csv_file_path = self.thread.csv_path.replace(".csv", "_data.csv")

        if os.path.exists(data_saving_csv_file_path):
            os.remove(data_saving_csv_file_path)
        
        self.prev_data_dict = {}

    def objective_function(self, x):
        """
        Evaluate the objective function for the optimization.
        
        This function runs the DATCOM simulation to compute a performance metric
        based on the aircraft geometry represented by x. It also updates the plot.
        If the evaluation exceeds a time limit or if a stop signal is detected, exceptions are raised.
        
        Args:
            x (list): Current values of the design variables.
        
        Returns:
            float: The computed wing load value, which serves as the objective function.
        
        Raises:
            TimeoutError: If the evaluation takes longer than 60 seconds.
            StopIteration: If a stop signal is detected before maximum iterations.
        """
        if (time.time() - self.t_iter) > 60:
            raise TimeoutError("Time limit exceeded while evaluating the objective function.")        
        elif self.thread_stop_event.is_set() and self.iteration != self.max_iter:
            self.iteration = self.max_iter
            raise StopIteration("Stopped by user")
        
        try:
            data, geometry = self.thread.loop(
                [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
                 x[3], x[4], x[5], x[6], x[7]],
                self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=None
            )
            data, geometry = self.thread.loop(
                [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
                 x[3], x[4], x[5], x[6], x[7]],
                self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=data["trim"]
            )
            self.prev_data_dict = data.copy()
            self.plot_class.plot_aircraft(geometry, x, self.iteration)

            E = self.methods.interpolate_fun(data, 'Alpha_list', 'Alpha_CL_list', self.alpha_cruise) / \
                self.methods.interpolate_fun(data, 'Alpha_list', 'Alpha_Cd_list', self.alpha_cruise)

            wing_load = self.mass * 9.81 / data["Sref"]

        except Exception as e:
            print("Error at evaluating the objective function: ", e)
            wing_load = 90
            E = 0
        return wing_load

    def constraints(self, x):
        """
        Evaluate the constraints for the optimization.
        
        This method runs DATCOM to obtain simulation data, computes stability derivatives
        and flight quality metrics, and assembles a list of inequality constraints.
        If the evaluation exceeds a time limit or if a stop signal is detected, exceptions are raised.
        
        Args:
            x (list): Current values of the design variables.
        
        Returns:
            list: A list of constraint values that must be non-negative.
        
        Raises:
            TimeoutError: If the evaluation takes longer than 60 seconds.
            StopIteration: If a stop signal is detected.
        """
        if (time.time() - self.t_iter) > 60:
            raise TimeoutError("Time limit exceeded while evaluating the constraint function.")
        elif self.thread_stop_event.is_set() and self.iteration != self.max_iter:
            raise StopIteration("Stopped by user")
        
        try:
            data, geometry = self.thread.loop(
                [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
                 x[3], x[4], x[5], x[6], x[7]],
                self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=None
            )
            iwb = data["trim"][0]
            it = data["trim"][1]
            data, geometry = self.thread.loop(
                [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
                 x[3], x[4], x[5], x[6], x[7]],
                self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=data["trim"]
            )
            self.prev_data_dict = data.copy()
            self.plot_class.plot_aircraft(geometry, x, self.iteration)

            derivatives = self.methods.get_stability_derivatives(data, geometry, self.alpha_cruise, self.flight_conditions[0], x)
            flight_quality = self.methods.get_flight_quality(data, derivatives, geometry, self.flight_conditions, self.mass, x[4])

            try:
                static_margin = -derivatives["Cma"] / derivatives["CLa"]
            except Exception:
                print("Error: CLa = 0")
                static_margin = -derivatives["Cma"]

            x_cg = x[0]
            fuselage_len = x[1]
            chrdr_w = x[2]
            chrdt_w = x[3]
            sspn_w = x[4]
            chrdr_htp = x[5]
            sspn_htp = x[6]
            x_w = x[7]

            basic_constraints = [
                -derivatives["Cma"],
                -derivatives["Cmq"],
                -derivatives["Clb"],
                -derivatives["Clp"],
                -derivatives["Cnr"],
                static_margin - 0.1,
                0.3 - static_margin,
            ]
            
            self.added_constraints = []
            for group_name, limits in self.default_limits.items():
                if group_name in ["Longitudinal", "Lateral-directional"]:
                    for group_name_mode, group_data in limits.items():
                        variables = group_data["variables"]
                        group_limits = group_data["limits"]

                        for idx, variable in enumerate(variables):
                            lower_limit = group_limits["lower"][idx]
                            upper_limit = group_limits["upper"][idx]

                            if group_name_mode == "Short Period":
                                mode = "CP"
                            elif group_name_mode == "Phugoid":
                                mode = "FUG"
                            elif group_name_mode == "Spiral":
                                mode = "SP"
                            elif group_name_mode == "Roll Mode":
                                mode = "CB"
                            elif group_name_mode == "Dutch Roll":
                                mode = "BH"

                            if lower_limit != "":
                                print("Lower", f"{mode}, {list(flight_quality[mode].keys())[idx]} > {float(lower_limit)} {flight_quality[mode][list(flight_quality[mode].keys())[idx]]}")
                                self.added_constraints.append(flight_quality[mode][list(flight_quality[mode].keys())[idx]] - float(lower_limit))

                            if upper_limit != "":
                                print("Higher", f"{mode}, {float(upper_limit)} > {list(flight_quality[mode].keys())[idx]} {flight_quality[mode][list(flight_quality[mode].keys())[idx]]}")
                                self.added_constraints.append(float(upper_limit) - flight_quality[mode][list(flight_quality[mode].keys())[idx]])

            constraints = basic_constraints + self.added_constraints
            E = self.methods.interpolate_fun(data, 'Alpha_list', 'Alpha_CL_list', self.alpha_cruise) / \
                self.methods.interpolate_fun(data, 'Alpha_list', 'Alpha_Cd_list', self.alpha_cruise)

            print("E = ", f"{E}"[:4], ", SM = ", f"{static_margin}"[:4], ", WL = ", f"{self.mass*9.81/data['Sref']}"[:4],
                  [constraint >= 0 for constraint in constraints])
            # Return flight quality constraints as a list of inequality constraints
            return constraints
        except Exception as e:
            print("Error at evaluating the constraint function: ", e)
            return -1

    def optimize(self):
        """
        Run the aerodynamic optimization loop.
        
        This method iteratively calls the objective function and constraint functions using SLSQP.
        It adjusts the bounds based on sensitivity until a solution is found or a termination condition is met.
        At the end, it prints final simulation data, derivatives, and flight quality metrics.
        """
        self.iteration = 0
        success_x = self.x0
        sensitivity = [0.05, 0.02, 0.01, 0.005]
        n = 0
        for k in range(self.max_iter):
            self.t_iter = time.time()

            # Adjust the bounds based on sensitivity
            new_bounds = []
            for i, (lower, upper) in enumerate(self.initial_bounds):
                lower_limit = max(self.x0[i] * (1 - sensitivity[n]), lower)
                upper_limit = min(self.x0[i] * (1 + sensitivity[n]), upper)
                new_bounds.append((lower_limit, upper_limit))
            self.bounds = new_bounds

            try:
                result = minimize(
                    self.objective_function, self.x0, method='SLSQP', bounds=self.bounds,
                    constraints={'type': 'ineq', 'fun': self.constraints},
                    options={'maxiter': self.max_iter,
                             'eps': self.eps if self.iteration != 0 else 1e-1,
                             'disp': True}
                )
                self.iteration += 1

                if result.success:
                    print("SUCCESS")
                    n = 0
                    self.x0 = result.x
                    success_x = result.x
                    continue
                else:
                    print("FINAL SOLUTION")
                    x = success_x
                    break

            except TimeoutError as e:
                print("TIME EXCEEDED")
                n += 1
                x = success_x
                if n == len(sensitivity):
                    print("FINAL SOLUTION")
                    break
            except StopIteration as e:
                print("STOPPED BY USER")
                x = success_x
                print(x)
                break

        self.iteration = self.max_iter
        data, geometry = self.thread.loop(
            [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
             x[3], x[4], x[5], x[6], x[7]],
            self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=None
        )
        data, geometry = self.thread.loop(
            [self.flight_conditions[0], self.flight_conditions[1], x[0], x[1], x[2],
             x[3], x[4], x[5], x[6], x[7]],
            self.flight_conditions, self.alpha_cruise, self.mass, self.prev_data_dict, to_trim=data["trim"]
        )
        self.plot_class.plot_aircraft(geometry, x, self.iteration)
        self.t_iter = time.time()
        derivatives = self.methods.get_stability_derivatives(data, geometry, self.alpha_cruise, self.flight_conditions[0], x)
        flight_quality = self.methods.get_flight_quality(data, derivatives, geometry, self.flight_conditions, self.mass, x[4])
        constraints = self.constraints(x)
        print("Constraints:", ", ".join([str(constraint) for constraint in constraints]))
        print("Objective function: ", self.objective_function(x))
        print("Final x: ", x)
        print(data)
        print(geometry)
        print(derivatives)
        print(flight_quality)
