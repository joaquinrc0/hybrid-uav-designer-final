"""
This module provides a Tkinter-based graphical interface for UAV redesign.
It allows users to input flight conditions, set variable limits, and run aerodynamic optimization.
Results are displayed both as textual updates and via a Matplotlib plot.
"""

import math
import threading
import os

import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from ambiance import Atmosphere

from scripts.plotter import ACPlotter
from scripts.optimizer import AerodynamicOptimization


class TkinterInterface:
    """
    Class that creates and manages the Tkinter interface for UAV redesign.
    
    Attributes:
        master (tk.Tk): The main Tkinter window.
        current_dir (str): The current working directory.
        style (ttk.Style): Custom style for widgets.
        notebook (ttk.Notebook): Notebook widget to hold different tabs.
        lower_limit (dict): Dictionary mapping variable names to their lower limit entry widgets.
        middle_value (dict): Dictionary mapping variable names to their functional value entry widgets.
        upper_limit (dict): Dictionary mapping variable names to their upper limit entry widgets.
        variable_labels (list): List of labels for displaying variable names.
        default_limits (dict): Dictionary containing default limits for variables.
        main_var_entries (list): List of entry widgets for common flight condition parameters.
        progress_bar (ttk.Progressbar): Progress bar widget for optimization progress.
        result_label (tk.Label): Label widget to display current status and results.
        variable_names (list): List of variable names used in result display.
        fig (matplotlib.figure.Figure): Figure for plotting optimization results.
        ax (matplotlib.axes.Axes): Axes for plotting optimization results.
        canvas (FigureCanvasTkAgg): Tkinter-compatible canvas for the Matplotlib figure.
    """

    def __init__(self, master):
        """
        Initialize the Tkinter interface, configure the main window and load or set default variable limits.

        Args:
            master (tk.Tk): The main Tkinter window.
        """
        self.master = master
        self.master.title("UAV Redesigner")
        self.master.state('zoomed')
        self.master.configure(bg='#1E1E1E')
        self.current_dir = os.getcwd()

        self.style = ttk.Style()
        self.style.configure('Custom.TButton', foreground='white', background='#4CAF50', font=('Helvetica', 12))
        self.style.configure('Custom.TLabelframe', background='#2E2E2E', font=('Helvetica', 12))
        self.style.configure('Custom.TNotebook', background='#2E2E2E')
        self.style.configure('Custom.TNotebook.Tab', background='#2E2E2E', foreground='black', font=('Helvetica', 14)) 

        self.notebook = ttk.Notebook(self.master, style='Custom.TNotebook')
        self.notebook.pack(fill="both", expand=True)

        self.lower_limit = {}
        self.middle_value = {}
        self.upper_limit = {}
        self.variable_labels = []

        # Define default limits from file or set a default configuration
        if os.path.exists(os.path.join(self.current_dir, "data", "last_config.txt")):
            with open(os.path.join(self.current_dir, "data", "last_config.txt"), "r") as file:
                self.default_limits = eval(file.readline())
        else:
            self.default_limits = {
                "General": {
                    "Geometry": {
                        "variables": [
                            "x_cg [m] (Mass center)",
                            "l [m] (Fuselage length)",
                            "chrdr_w [m] (Root chord of wing)",
                            "chrdt_w [m] (Tip chord of wing)",
                            "sspn_w [m] (Wing semi span)",
                            "chrdr_htp [m] (Root chord of horizontal tail plane)",
                            "sspn_htp [m] (Horizontal tail plane semi span)",
                            "x_w [m] (Wing position)"
                        ],
                        "limits": {
                            "lower": ["0.01", "0.1", "0.05", "0.05", "0.1", "0.05", "0.05", "0.1", "0.01", "0.01", "0.01"],
                            "middle": ["100", "100", "10", "10", "50", "10", "50", "100000", "100", "100", "100"],
                            "upper": ["100", "100", "10", "10", "50", "10", "50", "100000", "100", "100", "100"]
                        }
                    }
                },
                "Longitudinal": {
                    "Short Period": {
                        "variables": ["wn_SP", "damp_SP", "t1/2_SP"],
                        "limits": {
                            "lower": ["", "", "", "", "", ""],
                            "upper": ["", "", "", "", "", ""]
                        }
                    },
                    "Phugoid": {
                        "variables": ["damp_PH", "t1/2_PH"],
                        "limits": {
                            "lower": ["", "", "", "", "", ""],
                            "upper": ["", "", "", "", "", ""]
                        }
                    }
                },
                "Lateral-directional": {
                    "Spiral": {
                        "variables": ["t2_ME"],
                        "limits": {
                            "lower": [""],
                            "upper": [""]
                        }
                    },
                    "Roll Mode": {
                        "variables": ["tau_RM"],
                        "limits": {
                            "lower": [""],
                            "upper": [""]
                        }
                    },
                    "Dutch Roll": {
                        "variables": ["wn_BH", "damp_BH", "t1/2_BH"],
                        "limits": {
                            "lower": ["", "", "", ""],
                            "upper": ["", "", "", ""]
                        }
                    }
                }
            }

        # Create tabs for variable groups
        for group_name, limits in self.default_limits.items():
            self.create_tab(self.notebook, group_name, limits)
        
        self.create_result_tab(self.notebook, "Real Time Results")

        # Set up closing protocol for the main window
        self.master.protocol("WM_DELETE_WINDOW", self.close_window)

    def create_common_widgets(self, frame):
        """
        Create common input widgets for flight conditions.

        This method creates input fields for speed, height, cruise angle (alpha),
        maximum iterations, EPS value, and mass. It also adds buttons for running
        the optimization and stopping iterations.

        Args:
            frame (tk.Frame): The parent frame where the widgets are added.
        """
        row_frame = ttk.LabelFrame(frame, text="Flight Conditions", relief="groove", style='Custom.TLabelframe')
        row_frame.pack(padx=20, pady=5)

        self.main_var_entries = []

        # Speed input field
        tk.Label(row_frame, text="Speed [m/s]:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=0, column=0, padx=5, pady=5, sticky='w')
        self.speed_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.speed_entry.grid(row=0, column=1, padx=5, pady=5, sticky='e')
        self.speed_entry.insert(0, "43")
        self.main_var_entries.append(self.speed_entry)

        # Height input field
        tk.Label(row_frame, text="Height [m]:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=1, column=0, padx=5, pady=5, sticky='w')
        self.height_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.height_entry.grid(row=1, column=1, padx=5, pady=5, sticky='e')
        self.height_entry.insert(0, "1000")
        self.main_var_entries.append(self.height_entry)

        # Alpha Cruise input field
        tk.Label(row_frame, text="Alpha Cruise [°]:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=2, column=0, padx=5, pady=5, sticky='w')
        self.alpha_cruise_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.alpha_cruise_entry.grid(row=2, column=1, padx=5, pady=5, sticky='e')
        self.alpha_cruise_entry.insert(0, "2")
        self.main_var_entries.append(self.alpha_cruise_entry)

        # Iterations input field
        tk.Label(row_frame, text="Max Iterations:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=3, column=0, padx=5, pady=5, sticky='w')
        self.iterations_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.iterations_entry.grid(row=3, column=1, padx=5, pady=5, sticky='e')
        self.iterations_entry.insert(0, "20")
        self.main_var_entries.append(self.iterations_entry)

        # EPS input field
        tk.Label(row_frame, text="EPS:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=4, column=0, padx=5, pady=5, sticky='w')
        self.eps_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.eps_entry.grid(row=4, column=1, padx=5, pady=5, sticky='e')
        self.eps_entry.insert(0, "1e-2")
        self.main_var_entries.append(self.eps_entry)

        # Mass input field
        tk.Label(row_frame, text="Mass [kg]:", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
            .grid(row=5, column=0, padx=5, pady=5, sticky='w')
        self.mass_entry = tk.Entry(row_frame, font=('Helvetica', 14))
        self.mass_entry.grid(row=5, column=1, padx=5, pady=5, sticky='e')
        self.mass_entry.insert(0, "1.1")
        self.main_var_entries.append(self.mass_entry)

        # Create a frame for buttons
        button_frame = tk.Frame(frame)
        button_frame.pack(side=tk.BOTTOM, padx=5, pady=5)

        # Button to run results
        run_button = tk.Button(button_frame, text="Run Results", bg='#4CAF50', fg='white',
                               font=('Helvetica', 14), command=self.run_results)
        run_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Button to stop iterations
        stop_button = tk.Button(button_frame, text="Stop Iterations", bg='red', fg='white',
                                font=('Helvetica', 14), command=self.stop_iter)
        stop_button.pack(side=tk.LEFT, padx=5, pady=5)

    def stop_iter(self):
        """
        Signal the optimization thread to stop.

        This method sets the event flag used by the optimizer thread to signal a stop.
        """
        self.optimizer.thread_stop_event.set()

    def run_results(self):
        """
        Gather input values and launch the aerodynamic optimization.

        Reads variable limits and flight conditions from the input widgets, writes
        the configuration to file, and starts the optimization process in a new thread.
        """
        # Retrieve input values for "General" variables
        general_vars = [
            "x_cg [m] (Mass center)", "l [m] (Fuselage length)",
            "chrdr_w [m] (Root chord of wing)", "chrdt_w [m] (Tip chord of wing)",
            "sspn_w [m] (Wing semi span)", "chrdr_htp [m] (Root chord of horizontal tail plane)",
            "sspn_htp [m] (Horizontal tail plane semi span)", "x_w [m] (Wing position)"]
        
        general_limits = {"variables": general_vars, "limits": {"lower": [], "middle": [], "upper": []}}
        for var in general_vars:
            lower_limit = self.lower_limit[var].get()
            upper_limit = self.upper_limit[var].get()
            middle_value = self.middle_value[var].get()

            general_limits["limits"]["lower"].append(lower_limit)
            general_limits["limits"]["upper"].append(upper_limit)
            general_limits["limits"]["middle"].append(middle_value)

        # Retrieve input values for "Longitudinal" variables
        longitudinal_limits = {
            "Short Period": {
                "variables": ["wn_SP", "damp_SP", "t1/2_SP"],
                "limits": {"lower": [], "upper": []}
            },
            "Phugoid": {
                "variables": ["damp_PH", "t1/2_PH"],
                "limits": {"lower": [], "upper": []}
            }
        }
        for mode, values in longitudinal_limits.items():
            for var in values["variables"]:
                lower_limit = self.lower_limit[var].get()
                upper_limit = self.upper_limit[var].get()
                values["limits"]["lower"].append(lower_limit)
                values["limits"]["upper"].append(upper_limit)

        # Retrieve input values for "Lateral-directional" variables
        lateral_limits = {
            "Spiral": {
                "variables": ["t2_ME"],
                "limits": {"lower": [], "upper": []}
            },
            "Roll Mode": {
                "variables": ["tau_RM"],
                "limits": {"lower": [], "upper": []}
            },
            "Dutch Roll": {
                "variables": ["wn_BH", "damp_BH", "t1/2_BH"],
                "limits": {"lower": [], "upper": []}
            }
        }
        for mode, values in lateral_limits.items():
            for var in values["variables"]:
                lower_limit = self.lower_limit[var].get()
                upper_limit = self.upper_limit[var].get()
                values["limits"]["lower"].append(lower_limit)
                values["limits"]["upper"].append(upper_limit)

        # Update default limits and save configuration
        self.default_limits = {
            "General": {"Geometry": general_limits},
            "Longitudinal": longitudinal_limits,
            "Lateral-directional": lateral_limits
        }
        with open(os.path.join(self.current_dir, "data", "last_config.txt"), "w") as file:
            file.write(str(self.default_limits))

        # Prepare geometry dictionary and initialize the plotter
        geometry_dict = {
            'XCG_GEN': 0.6, 'ZCG_GEN': 0.0, 'XW_GEN': 0.5, 'ZW_GEN': 0.0,
            'ALIW_GEN': 8.17067, 'XH_GEN': 0.95, 'ZH_GEN': 0.0, 'ALIH_GEN': 12.1707,
            'XV_GEN': 0.97, 'ZV_GEN': 0.0, 'X(1)_BODY': [0.0, 0.009, 0.019, 0.033, 0.046, 0.054, 0.072, 0.113, 0.16, 0.654, 0.694, 0.791, 0.879, 0.943, 1.0],
            'NX_BODY': 15.0, 'ZU(1)_BODY': [0.0, 0.015, 0.021, 0.028, 0.04, 0.047, 0.056, 0.065, 0.069, 0.069, 0.069, 0.067, 0.064, 0.06, 0.05],
            'ZL(1)_BODY': [-0.0, -0.015, -0.02, -0.025, -0.027, -0.029, -0.031, -0.035, -0.036, -0.036, -0.034, -0.019, 0.0, 0.018, 0.036],
            'R(1)_BODY': [0.0, 0.015, 0.021, 0.027, 0.033, 0.036, 0.042, 0.048, 0.05, 0.05, 0.05, 0.045, 0.028, 0.02, 0.005],
            'CHRDR_WING': 0.2, 'CHRDTP_WING': 0.1, 'SSPN_WING': 0.3, 'SSPNE_WING': 0.27, 'CHSTAT_WING': 0.25,
            'TWISTA_WING': 0.0, 'TYPE_WING': 1.0, 'SAVSI_WING': 25.0,
            'CHRDR_HTP': 0.1, 'CHRDTP_HTP': 0.05, 'SSPN_HTP': 0.15, 'SSPNE_HTP': 0.135, 'SAVSI_HTP': 37.0,
            'CHSTAT_HTP': 0.0, 'TWISTA_HTP': 0.0, 'DHDADI_HTP': 0.0,
            'FTYPE_ELEVATOR': 1.0, 'NDELTA_ELEVATOR': 5.0,
            'DELTA(1)_ELEVATOR': [-20.0, -10.0, 0.0, 10.0, 20.0],
            'PHETE_ELEVATOR': 0.084, 'PHETEP_ELEVATOR': 0.084,
            'CHRDFI_ELEVATOR': 0.033, 'CHRDFO_ELEVATOR': 0.01,
            'SPANFI_ELEVATOR': 0.032, 'SPANFO_ELEVATOR': 0.148,
            'CB_ELEVATOR': 0.13048, 'TC_ELEVATOR': 0.05118,
            'CHRDR_VTP': 0.177, 'CHRDTP_VTP': 0.063, 'SSPN_VTP': 0.168, 'SSPNE_VTP': 0.168,
            'SAVSI_VTP': 40.0, 'CHSTAT_VTP': 0.0
        }
        self.plotter = ACPlotter(self.master, geometry_dict, self.variable_labels, self.variable_names,
                                   self.progress_bar, self.main_var_entries, self.fig, self.ax)

        # Set initial flight conditions
        velocity = float(self.main_var_entries[0].get())  # [m/s]
        alt = float(self.main_var_entries[1].get())         # [m]
        atm = Atmosphere(alt)

        # Prepare variable bounds and initial guess
        bounds = []
        x0 = []
        for var in self.default_limits["General"]["Geometry"]["variables"]:
            lower_limit = self.lower_limit[var].get()
            upper_limit = self.upper_limit[var].get()
            middle_value = self.middle_value[var].get()
            lower_limit = float(lower_limit) if lower_limit != "" else None
            upper_limit = float(upper_limit) if upper_limit != "" else None
            middle_value = float(middle_value) if middle_value != "" else None
            bounds.append((lower_limit, upper_limit))
            x0.append(middle_value)

        # Flight conditions tuple: (Mach, altitude, air density, velocity)
        flight_conditions = velocity / atm.speed_of_sound[0], alt, atm.density[0], velocity
        alpha_cruise = float(self.main_var_entries[2].get())
        
        self.thread_stop_event = threading.Event()
        self.thread = threading.Thread(target=self.optimize_helper, args=(
            self.thread_stop_event, x0, bounds, self.plotter, flight_conditions,
            alpha_cruise, self.main_var_entries,))
        self.thread.daemon = True
        self.thread.start()

    def create_tab(self, notebook, label_text, limits_dict):
        """
        Create a tab in the notebook for a specific group of variables.

        The tab contains labeled frames and entry widgets for setting lower and upper
        limits (and functional values for the General group).

        Args:
            notebook (ttk.Notebook): The notebook widget to add the tab.
            label_text (str): The title of the tab (e.g., "General", "Longitudinal").
            limits_dict (dict): Dictionary containing variable names and their default limits.
        """
        tab_frame = tk.Frame(notebook, bg='#1E1E1E')  # Tab background color
        tk.Label(tab_frame, text=label_text, bg='#1E1E1E', fg='white', font=('Helvetica', 16)).pack(pady=10)
        notebook.add(tab_frame, text=label_text)

        if label_text == "General":
            self.create_common_widgets(tab_frame)

        prev_group_name = ""
        for group_name, group_data in limits_dict.items():
            variables = group_data["variables"]
            group_limits = group_data["limits"]

            for idx, variable in enumerate(variables):
                if group_name != prev_group_name:
                    prev_group_name = group_name
                    row_frame = ttk.LabelFrame(tab_frame, text=group_name, relief="groove", style='Custom.TLabelframe')
                    row_frame.pack(padx=20, pady=5)

                variable_frame = tk.Frame(row_frame, bg='#1E1E1E')  # Background for variable frame
                variable_frame.pack(side=tk.TOP, padx=5, pady=5, fill="x", expand=True)
                for k in range(5): 
                    variable_frame.grid_columnconfigure(k, weight=1)

                self.lower_limit[variable] = tk.Entry(variable_frame, font=('Helvetica', 14))
                self.lower_limit[variable].insert(0, group_limits["lower"][idx])  # Insert default lower limit
                self.lower_limit[variable].grid(row=0, column=0, padx=5, pady=5, sticky="w")

                tk.Label(variable_frame, text=f"≤ {variable} ≤", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
                    .grid(row=0, column=1, padx=5, pady=5, sticky="nswe")
                
                self.upper_limit[variable] = tk.Entry(variable_frame, font=('Helvetica', 14))
                self.upper_limit[variable].insert(0, group_limits["upper"][idx])  # Insert default upper limit
                self.upper_limit[variable].grid(row=0, column=2, padx=5, pady=5, sticky="e")

                if label_text == "General":
                    tk.Label(variable_frame, text=f"Functional value →", bg='#1E1E1E', fg='white', font=('Helvetica', 14))\
                        .grid(row=0, column=3, pady=5, sticky="e")
                    self.middle_value[variable] = tk.Entry(variable_frame, font=('Helvetica', 14))
                    self.middle_value[variable].insert(0, group_limits["middle"][idx])  # Insert default middle value
                    self.middle_value[variable].grid(row=0, column=4, padx=5, pady=5, sticky="e")

    def create_result_tab(self, notebook, label_text):
        """
        Create the result tab for displaying optimization outputs and graphs.

        This tab features a progress bar, status label, variable labels, and a Matplotlib
        plot embedded in a Tkinter-compatible canvas.

        Args:
            notebook (ttk.Notebook): The notebook widget to add the tab.
            label_text (str): The title of the results tab.
        """
        result_tab_frame = tk.Frame(notebook, bg='#1E1E1E')
        notebook.add(result_tab_frame, text=label_text)

        for i in range(4):
            result_tab_frame.grid_columnconfigure(i, weight=1)
        result_tab_frame.grid_rowconfigure(0, weight=1)

        # Create subframes for layout
        left_frame = tk.Frame(result_tab_frame, bg='#1E1E1E')
        left_frame.grid(row=0, column=0, padx=(10, 50), pady=10, sticky="nswe")
        left_frame.grid_columnconfigure(0, weight=1)
        left_frame.grid_rowconfigure(0, weight=1)

        self.right_frame = tk.Frame(result_tab_frame, bg='#1E1E1E')
        self.right_frame.grid(row=0, column=1, columnspan=3, padx=(10, 50), pady=10, sticky="nswe")
        self.right_frame.grid_columnconfigure(0, weight=1)
        self.right_frame.grid_rowconfigure(0, weight=1)

        # Progress bar for optimization progress
        self.progress_bar = ttk.Progressbar(left_frame, orient='horizontal', length=200, mode='determinate')
        self.progress_bar.grid(row=0, column=0, padx=10, pady=10, sticky="we")

        # Label for status updates
        self.result_label = tk.Label(left_frame, text="", bg='#1E1E1E', fg='white', font=('Helvetica', 14))
        self.result_label.grid(row=1, column=0, padx=10, pady=10, sticky=tk.N)

        # Variable labels for result display
        self.variable_names = ["STATUS", "X_CG", "FUSELAGE_LEN", "CHRDR_W", "CHRDT_W",
                               "SSPN_W", "CHRDR_HTP", "SSPN_HTP", "X_W"]
        self.variable_labels = []

        for idx, variable in enumerate(self.variable_names):
            label = tk.Label(left_frame, text=variable, bg='#1E1E1E', fg='white',
                             font=('Helvetica', 14, "bold"))
            label.grid(row=idx + 1, column=0, padx=10, pady=5, sticky="e")
            self.variable_labels.append(label)

        # Configure row weights
        for i in range(len(self.variable_names) + 1):
            left_frame.grid_rowconfigure(i, weight=1)

        # Matplotlib figure for plotting results
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky="nswe")

        left_frame.grid_propagate(False)
        self.right_frame.grid_propagate(False)

    def optimize_helper(self, thread_stop_event, x0, bounds, plotter, flight_conditions, alpha_cruise, main_var_entries):
        """
        Helper method to run aerodynamic optimization in a separate thread.

        This method creates an instance of AerodynamicOptimization with the provided parameters
        and starts the optimization process.

        Args:
            thread_stop_event (threading.Event): Event to signal thread termination.
            x0 (list): Initial guess for optimization variables.
            bounds (list): List of tuples specifying lower and upper bounds for each variable.
            plotter (ACPlotter): Plotter object for updating the GUI with results.
            flight_conditions (tuple): Flight condition parameters (Mach, altitude, density, velocity).
            alpha_cruise (float): Cruise angle of attack.
            main_var_entries (list): List of Tkinter entry widgets for flight parameters.
        """
        self.optimizer = AerodynamicOptimization(
            thread_stop_event, x0, bounds, plotter, flight_conditions, alpha_cruise, main_var_entries, self.default_limits)
        self.optimizer.optimize()

    def close_window(self):
        """
        Handle closing of the Tkinter interface.

        Closes the Matplotlib figure and destroys the main Tkinter window.
        """
        plt.close()  # Close the Matplotlib figure
        self.master.destroy()  # Destroy the main Tkinter window


def main():
    """
    Main function to launch the UAV Redesigner interface.

    Creates the root Tkinter window, instantiates the TkinterInterface, and starts the main event loop.
    """
    root = tk.Tk()
    app = TkinterInterface(root)
    root.mainloop()


if __name__ == "__main__":
    main()
