import matplotlib.pyplot as plt
import math

class ACPlotter:
    """
    A class for plotting the aircraft geometry during aerodynamic optimization.

    This class uses Matplotlib to visualize various components of the aircraft,
    such as the fuselage, wing, horizontal and vertical stabilizers, and center of gravity.
    It also updates a progress bar and variable labels based on the current optimization iteration.

    Attributes:
        master (tk.Widget): The parent Tkinter widget.
        geometry_dict (dict): Dictionary containing the aircraft geometry parameters.
        variable_labels (list): List of Tkinter label widgets to display variable information.
        variable_names (list): List of variable names corresponding to the labels.
        progress_bar (ttk.Progressbar): Progress bar widget to indicate optimization progress.
        max_iterations (int): Maximum number of iterations for the optimization process.
        fig (matplotlib.figure.Figure): Matplotlib figure object for plotting.
        ax (matplotlib.axes.Axes): Matplotlib axes object for plotting.
    """

    def __init__(self, master, geometry_dict, variable_labels, variable_names, progress_bar, main_var_entries, fig, ax):
        """
        Initialize the ACPlotter instance.

        Args:
            master (tk.Widget): The parent Tkinter widget.
            geometry_dict (dict): Initial dictionary with aircraft geometry parameters.
            variable_labels (list): List of Tkinter label widgets for displaying variable data.
            variable_names (list): List of variable names used for labeling.
            progress_bar (ttk.Progressbar): Progress bar widget to show optimization progress.
            main_var_entries (list): List of entry widgets containing main parameters; used here to extract the maximum iterations.
            fig (matplotlib.figure.Figure): Matplotlib figure for plotting.
            ax (matplotlib.axes.Axes): Matplotlib axes for plotting.
        """
        self.master = master
        self.geometry_dict = geometry_dict
        self.variable_labels = variable_labels
        self.variable_names = variable_names
        self.progress_bar = progress_bar
        self.max_iterations = int(main_var_entries[3].get())
        self.fig, self.ax = fig, ax

    def plot_aircraft(self, geometry_dict=None, x=None, iteration=None):
        """
        Plot the aircraft geometry and update the progress bar and variable labels.

        If a new geometry_dict is provided, it updates the internal state and then plots the aircraft's
        fuselage, wing, horizontal stabilizer, vertical stabilizer, and center of gravity.

        Args:
            geometry_dict (dict, optional): New aircraft geometry parameters. If provided, replaces the current geometry.
            x (list, optional): Current design variable values.
            iteration (int, optional): The current iteration count for optimization.
        """
        if geometry_dict:
            self.geometry_dict = geometry_dict
            self.x = x
            max_len = 20
            self.progress_bar["value"] = int((iteration / self.max_iterations) * 100)
            for idx, var in enumerate(self.variable_labels):
                if idx != 0:
                    var["text"] = f"{self.variable_names[idx]: <{max_len}} = {str(x[idx-1])[:4]}"
                else:
                    if iteration == self.max_iterations:
                        var["text"] = f"{self.variable_names[idx]: <{max_len}} = STANDING"
                        var.config(fg="red")
                    else:
                        var["text"] = f"{self.variable_names[idx]: <{max_len}} = RUNNING..."
                        var.config(fg="green")

            self.ax.clear()

            # Extract general geometry variables
            XCG_GEN = geometry_dict["XCG_GEN"]
            ZCG_GEN = geometry_dict["ZCG_GEN"]
            XW_GEN = geometry_dict["XW_GEN"]
            ZW_GEN = geometry_dict["ZW_GEN"]
            ALIW_GEN = geometry_dict["ALIW_GEN"]
            XH_GEN = geometry_dict["XH_GEN"]
            ZH_GEN = geometry_dict["ZH_GEN"]
            ALIH_GEN = geometry_dict["ALIH_GEN"]
            XV_GEN = geometry_dict["XV_GEN"]
            ZV_GEN = geometry_dict["ZV_GEN"]

            # Extract body geometry variables
            X_BODY = geometry_dict["X(1)_BODY"]
            NX_BODY = geometry_dict["NX_BODY"]
            ZU_BODY = geometry_dict["ZU(1)_BODY"]
            ZL_BODY = geometry_dict["ZL(1)_BODY"]
            R_BODY = geometry_dict["R(1)_BODY"]

            # Extract wing geometry variables
            CHRDR_WING = geometry_dict["CHRDR_WING"]
            CHRDTP_WING = geometry_dict["CHRDTP_WING"]
            SSPN_WING = geometry_dict["SSPN_WING"]
            SSPNE_WING = geometry_dict["SSPNE_WING"]
            CHSTAT_WING = geometry_dict["CHSTAT_WING"]
            TWISTA_WING = geometry_dict["TWISTA_WING"]
            TYPE_WING = geometry_dict["TYPE_WING"]
            SAVSI_WING = geometry_dict["SAVSI_WING"]

            # Extract horizontal tail geometry variables
            CHRDR_HTP = geometry_dict["CHRDR_HTP"]
            CHRDTP_HTP = geometry_dict["CHRDTP_HTP"]
            SSPN_HTP = geometry_dict["SSPN_HTP"]
            SSPNE_HTP = geometry_dict["SSPNE_HTP"]
            SAVSI_HTP = geometry_dict["SAVSI_HTP"]
            CHSTAT_HTP = geometry_dict["CHSTAT_HTP"]
            TWISTA_HTP = geometry_dict["TWISTA_HTP"]
            DHDADI_HTP = geometry_dict["DHDADI_HTP"]

            # Extract elevator geometry variables
            FTYPE_ELEVATOR = geometry_dict["FTYPE_ELEVATOR"]
            NDELTA_ELEVATOR = geometry_dict["NDELTA_ELEVATOR"]
            DELTA_ELEVATOR = geometry_dict["DELTA(1)_ELEVATOR"]
            PHETE_ELEVATOR = geometry_dict["PHETE_ELEVATOR"]
            PHETEP_ELEVATOR = geometry_dict["PHETEP_ELEVATOR"]
            CHRDFI_ELEVATOR = geometry_dict["CHRDFI_ELEVATOR"]
            CHRDFO_ELEVATOR = geometry_dict["CHRDFO_ELEVATOR"]
            SPANFI_ELEVATOR = geometry_dict["SPANFI_ELEVATOR"]
            SPANFO_ELEVATOR = geometry_dict["SPANFO_ELEVATOR"]
            CB_ELEVATOR = geometry_dict["CB_ELEVATOR"]
            TC_ELEVATOR = geometry_dict["TC_ELEVATOR"]

            # Extract vertical tail geometry variables
            CHRDR_VTP = geometry_dict["CHRDR_VTP"]
            CHRDTP_VTP = geometry_dict["CHRDTP_VTP"]
            SSPN_VTP = geometry_dict["SSPN_VTP"]
            SSPNE_VTP = geometry_dict["SSPNE_VTP"]
            SAVSI_VTP = geometry_dict["SAVSI_VTP"]
            CHSTAT_VTP = geometry_dict["CHSTAT_VTP"]

            # Calculate fuselage coordinates
            x_fuselaje = [x for x in R_BODY]
            y_fuselaje = [-y for y in X_BODY]
            x_fuselaje.extend([-x for x in x_fuselaje])
            y_fuselaje.extend([y for y in y_fuselaje])

            # Calculate wing coordinates
            wing_x = [0, SSPN_WING, SSPN_WING, 0]
            wing_y = [
                -XW_GEN,
                -XW_GEN - SSPN_WING * math.tan(math.radians(SAVSI_WING)),
                -XW_GEN - SSPN_WING * math.tan(math.radians(SAVSI_WING)) - CHRDTP_WING,
                -XW_GEN - CHRDR_WING
            ]
            wing_x.extend([-x for x in reversed(wing_x)])
            wing_y.extend([y for y in wing_y])

            # Center of Gravity coordinates
            cg_x = [0]
            cg_y = [-XCG_GEN]  # Arbitrary CG y-position

            # Calculate horizontal stabilizer coordinates
            htp_x = [0, SSPN_HTP, SSPN_HTP, 0]
            htp_y = [
                -XH_GEN,
                -XH_GEN - SSPN_HTP * math.tan(math.radians(SAVSI_HTP)),
                -XH_GEN - SSPN_HTP * math.tan(math.radians(SAVSI_HTP)) - CHRDTP_HTP,
                -XH_GEN - CHRDR_HTP
            ]
            htp_x.extend([-x for x in reversed(htp_x)])
            htp_y.extend([y for y in htp_y])

            # Calculate vertical stabilizer coordinates
            vtp_x = [0, SSPN_VTP, SSPN_VTP, 0]
            vtp_y = [
                -XV_GEN,
                -XV_GEN - SSPN_VTP * math.tan(math.radians(SAVSI_VTP)),
                -XV_GEN - SSPN_VTP * math.tan(math.radians(SAVSI_VTP)) - CHRDTP_VTP,
                -XV_GEN - CHRDR_VTP
            ]

            # Plot the aircraft components
            self.ax.plot(wing_x, wing_y, 'r-', label='Wing')
            self.ax.plot(cg_x, cg_y, 'kx', markersize=10, label='Center of Gravity (CG)')
            self.ax.plot(htp_x, htp_y, 'g-', label='Horizontal Stabilizer')
            self.ax.plot(vtp_x, vtp_y, 'y-', label='Vertical Stabilizer')
            self.ax.plot(x_fuselaje, y_fuselaje, color='blue', label='Fuselage')
            self.ax.set_aspect('equal', adjustable='box')
            self.ax.set_xlabel('Width (m)')
            self.ax.set_ylabel('Length (m)')
            self.ax.set_title('Top View of the Aircraft')
            self.ax.legend(fontsize='small')
            self.ax.grid(True)

            plt.tight_layout()
            plt.draw()
            self.master.update()
