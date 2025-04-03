class DCMCreator:
    """
    A class to create and write a DCM (Digital Configuration Model) file for an aircraft configuration.

    This class calculates and formats various geometric and aerodynamic parameters of the aircraft
    based on the input values provided (such as fuselage length, wing geometry, etc.), and writes
    the configuration to a file in a specific format.

    :param filename: The output file name where the DCM configuration will be written.
    :param mach: Mach number for the flight condition.
    :param alt: Altitude for the flight condition.
    :param x_cg: x-coordinate of the aircraft's center of gravity.
    :param fuselage_len: Length of the fuselage, used for scaling geometric parameters.
    :param chrdr_w: Chord (root) dimension for the wing.
    :param chrdt_w: Chord (tip) dimension for the wing.
    :param sspn_w: Wing span used in the wing parameter definitions.
    :param chrdr_htp: Chord dimension for the horizontal tail plane.
    :param sspn_htp: Span of the horizontal tail plane.
    :param x_w: x-coordinate position of the wing.
    :param to_trim: A list of trim angles; the first element is for the wing and the second for the horizontal tail.
    """

    def __init__(self, filename, mach, alt, x_cg, fuselage_len, chrdr_w, chrdt_w, sspn_w, chrdr_htp, sspn_htp, x_w, to_trim):
        x_htp = 0.95 * fuselage_len
        x_vtp = 0.97 * fuselage_len
        self.filename = filename
        self.geometry_dict = {}
        if to_trim is None:
            to_trim = [0.0, 0.0]

        # Commands and parameters to be written in the DCM file
        self.commands_and_params = [
            "DIM M",
            "PART",
            "DAMP",
            "BUILD",
            "TRIM",
            "DERIV RAD",
        ]

        # Flight condition parameters
        self.fltcon_params = {
            "NMACH": 1.0,
            "MACH(1)": "{:.{}f}".format(mach, 3),
            "NALT": 1.0,
            "ALT(1)": "{:.{}f}".format(alt, 3),
            "NALPHA": 13.0,
            "ALSCHD(1)": [-5.0, 0.0, 1.9, 2.0, 2.1, 2.81, 3.0, 3.69, 4.5, 6.0, 8.0, 11.7, 13.0]
        }

        self.close_fltcon = f"  STMACH={0.8}$"

        # General synthesis parameters for aircraft geometry
        self.synths_params = {
            "XCG": x_cg,
            "ZCG": 0.0,
            "XW": x_w,
            "ZW": 0.0,
            "ALIW": to_trim[0],
            "XH": x_htp,
            "ZH": 0.0,
            "ALIH": to_trim[1],
            "XV": x_vtp,
            "ZV": 0.0,
        }

        self.close_synths = "  VERTUP=.TRUE.$"

        # Fuselage parameters and scaling
        fuselage_x_list = [0.000000, 0.349000, 0.750000, 1.312000, 1.803000, 2.124000,
                           2.839000, 4.479000, 6.350000, 25.889000, 27.459000, 31.289000,
                           34.789000, 37.296000, 39.570000]

        fuselage_zu_list = [0.010000, 0.590000, 0.850000, 1.110000, 1.590000, 1.860000,
                            2.230000, 2.580000, 2.730000, 2.730000, 2.720000, 2.660000, 2.520000,
                            2.390000, 1.959000]
        fuselage_zl_list = [-0.010000, -0.590000, -0.780000, -0.970000, -1.080000, -1.140000,
                            -1.240000, -1.370000, -1.420000, -1.420000, -1.340000, -0.770000,
                            0.000000, 0.703000, 1.422000]
        fuselage_r_list = [0.010000, 0.587000, 0.816000, 1.086000, 1.298000, 1.414000,
                           1.647000, 1.895000, 1.976000, 1.976000, 1.971000, 1.799000, 1.113000,
                           0.806000, 0.204000]
        self.body_params = {
            "X(1)": ["{:.{}f}".format(element * fuselage_len / 39.57, 3) for element in fuselage_x_list],
            "NX": 15.000000,
            "ZU(1)": ["{:.{}f}".format(element * fuselage_len / 39.57, 3) for element in fuselage_zu_list],
            "ZL(1)": ["{:.{}f}".format(element * fuselage_len / 39.57, 3) for element in fuselage_zl_list],
            "R(1)": ["{:.{}f}".format(element * fuselage_len / 39.57, 3) for element in fuselage_r_list],
        }

        self.close_body = "  METHOD=2.000000$"

        # Wing parameters
        self.wing_params = {
            "CHRDR": "{:.{}f}".format(chrdr_w, 3),
            "CHRDTP": "{:.{}f}".format(chrdt_w, 3),
            "SSPN": "{:.{}f}".format(sspn_w, 3),
            "SSPNE": "{:.{}f}".format(sspn_w - 1.97 * fuselage_len / 39.57, 3),
            "CHSTAT": "0.25",
            "TWISTA": 0.0,
            "TYPE": 1.0,
            "SAVSI": "{:.{}f}".format(8.34, 3)
        }

        self.close_wing = f"  DHDADI={5.0}$"
        self.naca_wing = "NACA-W-6-65-108A"

        # Ailerons parameters (commented out in the write method)
        self.ailerons_params = {
            "STYPE": 4.0,
            "NDELTA": 9.0,
            "DELTAL": [-20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0],
            "DELTAR": [20.0, 15.0, 10.0, 5.0, 0.0, -5.0, -10.0, -15.0, -20.0],
            "SPANFI": "{:.{}f}".format(11.935 * fuselage_len / 39.57, 3),
            "SPANFO": "{:.{}f}".format(16.1975 * fuselage_len / 39.5, 3),
            "CHRDFI": "{:.{}f}".format(0.727175 * fuselage_len / 39.57, 3)
        }

        self.close_ailerons = "  CHRDFO={:.{}f}".format(0.4422 * fuselage_len / 39.57, 3)

        # Horizontal Tail Plane (HTP) parameters
        self.htp_params = {
            "CHRDR": "{:.{}f}".format(chrdr_htp, 3),
            "CHRDTP": "{:.{}f}".format(chrdr_htp, 3),
            "SSPN": "{:.{}f}".format(sspn_htp, 3),
            "SSPNE": "{:.{}f}".format(sspn_htp - 0.8 * fuselage_len / 39.57, 3),
            "SAVSI": "{:.{}f}".format(37.0, 3),
            "CHSTAT": 0.0,
            "TWISTA": 0.0,
            "DHDADI": 0.0
        }

        self.close_htp = f"  TYPE={1.0}$"
        self.naca_htp = "NACA-H-6-66-008A"

        # Elevator parameters
        self.elevator_params = {
            "FTYPE": 1.0,
            "NDELTA": 5.0,
            "DELTA(1)": [-20.0, -10.0, 0.0, 10.0, 20.0],
            "PHETE": 0.084,
            "PHETEP": 0.084,
            "CHRDFI": "{:.{}f}".format(1.304775 * fuselage_len / 39.57, 3),
            "CHRDFO": "{:.{}f}".format(0.389075 * fuselage_len / 39.57, 3),
            "SPANFI": "{:.{}f}".format(1.254815 * fuselage_len / 39.57, 3),
            "SPANFO": "{:.{}f}".format(5.841485 * fuselage_len / 39.57, 3),
            "CB": 0.13048,
            "TC": 0.05118
        }

        self.close_elevator = f"  NTYPE={1.0}$"

        # Vertical Tail Plane (VTP) parameters
        self.vtp_params = {
            "CHRDR": "{:.{}f}".format(0.2 * fuselage_len, 3),
            "CHRDTP": "{:.{}f}".format(2.5 / 7.0 * 0.2 * fuselage_len, 3),
            "SSPN": "{:.{}f}".format(0.2 * fuselage_len, 3),
            "SSPNE": "{:.{}f}".format(0.2 * fuselage_len, 3),
            "SAVSI": 40.0,
            "CHSTAT": 0.0
        }

        self.close_vtp = f"  TYPE={1.0}$"
        self.naca_vtp = "NACA-V-6-66-008A"

    def write_dcm_file(self):
        """
        Writes the complete DCM configuration to the file specified by self.filename.

        This method writes different sections of the DCM file including flight condition,
        general synthesis, fuselage (body), wing, horizontal tail plane, elevator, and vertical tail plane.
        Each section is formatted appropriately with its corresponding parameters.

        :raises IOError: If the file cannot be opened or written.
        """
        with open(self.filename, 'w') as f:
            # Commands and parameters
            for command in self.commands_and_params:
                f.write(command + "\n")

            # Flight condition configuration
            f.write(" $FLTCON\n")
            for key, value in self.fltcon_params.items():
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_fltcon + "\n")

            # General synthesis configuration
            f.write(" $SYNTHS\n")
            for key, value in self.synths_params.items():
                self.geometry_dict[f"{key}_GEN"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_synths + "\n")

            # Fuselage (body) configuration
            f.write(" $BODY\n")
            for key, value in self.body_params.items():
                self.geometry_dict[f"{key}_BODY"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_body + "\n")

            # Wing configuration
            f.write(" $WGPLNF\n")
            for key, value in self.wing_params.items():
                self.geometry_dict[f"{key}_WING"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_wing + "\n")
            f.write(self.naca_wing + "\n")

            # Ailerons section is commented out
            # f.write(" $ASYFLP\n")
            # for key, value in self.ailerons_params.items():
            #     if isinstance(value, list):
            #         f.write(self.split_list(key, value))
            #     else:
            #         f.write("  " + key + "=" + str(value) + ",\n")
            # f.write(self.close_ailerons + "\n")

            # Horizontal Tail Plane (HTP) configuration
            f.write(" $HTPLNF\n")
            for key, value in self.htp_params.items():
                self.geometry_dict[f"{key}_HTP"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_htp + "\n")
            f.write(self.naca_htp + "\n")

            # Elevator configuration
            f.write(" $SYMFLP\n")
            for key, value in self.elevator_params.items():
                self.geometry_dict[f"{key}_ELEVATOR"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_elevator + "\n")

            # Vertical Tail Plane (VTP) configuration
            f.write(" $VTPLNF\n")
            for key, value in self.vtp_params.items():
                self.geometry_dict[f"{key}_VTP"] = value
                if isinstance(value, list):
                    f.write(self.split_list(key, value))
                else:
                    f.write("  " + key + "=" + str(value) + ",\n")
            f.write(self.close_vtp + "\n")
            f.write(self.naca_vtp + "\n")

    def split_list(self, key, input_list):
        """
        Formats a list of numerical values for inclusion in the DCM file.

        The input list is converted into a string of comma-separated values,
        grouped in chunks of three per line with proper indentation. This helps
        maintain readability in the DCM file.

        :param key: The parameter name associated with the list.
        :param input_list: A list of numbers to be formatted.
        :return: A formatted string representing the list, ready to be written to the DCM file.
        """
        input_list = [str(number) for number in input_list]
        result = []
        for i in range(0, len(input_list), 3):
            result.append("  " + ','.join(input_list[i:i+3]) + ",")
        return "  " + key + "=" + '\n'.join(result)[2:] + "\n"
