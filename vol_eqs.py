import pandas as pd
import numpy as np
from ctypes import *

# A dictionary to map error codes from the DLL to human-readable messages
# This mimics the 'GetError' and 'ErrorMessages' functions in the VBA script.
VOLLIB_ERROR_CODES = {
    1: "Invalid or missing Volume Equation string.",
    2: "Missing Form Class for a Form Class equation.",
    3: "DBH is less than 1.0 inch.",
    4: "Tree height is less than 4.5 feet.",
    5: "The D-squared-H value is out of range for the equation.",
    6: "Invalid species code for the equation.",
    7: "Illegal primary product log height.",
    8: "Illegal secondary product log height.",
    9: "Upper stem diameter measurement is required.",
    10: "Illegal upper stem height (less than 4.5ft or greater than total height).",
    11: "Unable to fit profile model with the given DBH, merch height, and top DIB.",
    12: "Tree has more than 20 logs.",
    13: "Top DIB is greater than DBH inside bark.",
    14: "Bark equation does not exist or yields negative DBHIB.",
    15: "Invalid Biomass Equation specified.",
    16: "Primary product height (HT1PRD) is required for this biomass calculation.",
    17: "Secondary product height (HT2PRD) is required for this biomass calculation.",
    18: "Main stem cubic volume (CV4) is needed for this biomass calculation.",
    19: "The specified equation is not a profile model.",
    20: "Invalid Region number.",
    21: "Invalid Forest number.",
    22: "Invalid District number.",
    23: "Invalid FIA species code.",
}

class VollibWrapper:
    """
    A Python wrapper for the NVEL (VOLLIB) FORTRAN library that mimics the
    functionality of the provided VBA script.

    This class loads the vollib.dll, defines the necessary function prototypes,
    and provides high-level methods to calculate various tree volumes, biomass
    weights, and taper information.

    Attributes:
        dll_path (str): The file path to the vollib.dll file.
        vollib (WinDLL): The loaded DLL object from ctypes.
    """
    def __init__(self, dll_path: str):
        """
        Initializes the wrapper by loading the DLL and setting up function signatures.

        Args:
            dll_path (str): The full path to the 'vollib.dll' or 'vollib64.dll' file.
        
        Raises:
            OSError: If the DLL cannot be loaded from the specified path.
        """
        self.dll_path = dll_path
        try:
            self.vollib = windll.LoadLibrary(self.dll_path)
            print(f"Successfully loaded vollib.dll from: {self.dll_path}")
        except OSError as e:
            print(f"Error loading DLL: {e}")
            print("Please ensure the path is correct and it matches your Python architecture (32-bit vs 64-bit).")
            raise

        self._define_function_prototypes()

    def _define_function_prototypes(self):
        """Sets the argument types (argtypes) for the DLL functions."""
        # For GETVOLEQ
        self.vollib.GETVOLEQ.argtypes = [
            POINTER(c_int), c_char_p, c_int, c_char_p, c_int, POINTER(c_int),
            c_char_p, c_int, c_char_p, c_int, POINTER(c_int)
        ]
        # For GETNVBEQ2 (used for biomass)
        self.vollib.GETNVBEQ2.argtypes = [
            POINTER(c_int), c_char_p, c_int, c_char_p, c_int, POINTER(c_int),
            c_char_p, c_int, POINTER(c_int)
        ]
        # For VOLUMELIBRARY2
        self.vollib.VOLUMELIBRARY2.argtypes = [
            POINTER(c_int), c_char_p, c_int, c_char_p, c_int, POINTER(c_float),
            POINTER(c_float), POINTER(c_float), POINTER(c_float), POINTER(c_float),
            c_char_p, c_int, POINTER(c_float), POINTER(c_int), POINTER(c_float),
            POINTER(c_float), POINTER(c_float), POINTER(c_float), POINTER(c_float),
            POINTER(c_float), POINTER(c_int), POINTER(c_float), POINTER(c_float),
            POINTER(c_int), POINTER(c_float), POINTER(c_float),
            # Output Arrays
            POINTER((c_float * 15)), POINTER((c_float * 7 * 20)),
            POINTER((c_float * 21 * 3)), POINTER((c_float * 20)),
            POINTER((c_float * 21)), POINTER(c_int), POINTER(c_float),
            POINTER(c_float), POINTER(c_int), POINTER(c_int), POINTER(c_int),
            POINTER(c_int), POINTER(c_int), c_char_p, c_int, c_char_p, c_int,
            POINTER(c_int), c_char_p, c_int, POINTER(c_int), POINTER(c_int),
            c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_float),
            POINTER(c_float), POINTER(c_int),
            # Biomass Output Arrays
            POINTER((c_float * 15)), POINTER((c_float * 15)),
            POINTER(c_float), POINTER(c_float), POINTER(c_int)
        ]

    def _get_volume_equation(self, region: int, forest: str, district: str,
                             species_code: int, eq_type: str = 'volume') -> str:
        """
        Internal method to retrieve the appropriate volume equation string.

        Args:
            region (int): Forest Service region number.
            forest (str): Forest code (2-digit string).
            district (str): District code (2-digit string).
            species_code (int): FIA species code.
            eq_type (str): 'volume' for standard volume, 'biomass' for NSVB.

        Returns:
            The volume equation string.
        
        Raises:
            ValueError: If the DLL call fails.
        """
        c_regn = c_int(region)
        c_forst = create_string_buffer(f'{forest:<2}'.encode('ascii'))
        c_dist = create_string_buffer(f'{district:<2}'.encode('ascii'))
        c_spec = c_int(species_code)
        c_voleq = create_string_buffer(11) # 10 chars + null terminator
        c_errflag = c_int(0)

        if eq_type == 'volume':
            c_prod = create_string_buffer(b'01')
            self.vollib.GETVOLEQ(
                byref(c_regn), c_forst, len(c_forst)-1, c_dist, len(c_dist)-1,
                byref(c_spec), c_prod, len(c_prod)-1, c_voleq, len(c_voleq)-1,
                byref(c_errflag)
            )
        elif eq_type == 'biomass':
             self.vollib.GETNVBEQ2(
                byref(c_regn), c_forst, len(c_forst)-1, c_dist, len(c_dist)-1,
                byref(c_spec), c_voleq, len(c_voleq)-1, byref(c_errflag)
            )
        else:
            raise ValueError("eq_type must be 'volume' or 'biomass'")

        if c_errflag.value != 0:
            error_message = VOLLIB_ERROR_CODES.get(c_errflag.value, f"Unknown error code {c_errflag.value} from GETVOLEQ/GETNVBEQ2.")
            raise ValueError(error_message)

        return c_voleq.value.decode('ascii').strip()

    def _call_volume_library(self, **kwargs):
        """
        The core engine method that prepares variables and calls VOLUMELIBRARY2.
        This is the central point for all public-facing calculation methods.
        """
        # --- 1. Set required and optional parameters from kwargs ---
        region = kwargs.get('region', 1)
        forest = kwargs.get('forest', '01')
        district = kwargs.get('district', '01')
        species_code = kwargs.get('species_code')
        dbh = kwargs.get('dbh')
        height = kwargs.get('height')

        voleq_override = kwargs.get('voleq_override')
        merch_top_dib_primary = kwargs.get('merch_top_dib_primary', 4.0)
        merch_top_dib_secondary = kwargs.get('merch_top_dib_secondary', 0.0)
        stump_height = kwargs.get('stump_height', 1.0)
        
        # FVS / Biomass specific parameters
        crown_ratio = kwargs.get('crown_ratio', 0.0)
        cull = kwargs.get('cull', 0.0)
        decay_cd = kwargs.get('decay_cd', 0)
        broken_ht = kwargs.get('broken_ht', 0.0)
        live_dead = kwargs.get('live_dead', 'L')
        calc_type = kwargs.get('calc_type', 'FVS')

        # --- 2. Get Volume Equation if not provided ---
        eq_type = kwargs.get('eq_type', 'volume')
        if voleq_override:
            voleq_str = voleq_override
        else:
            voleq_str = self._get_volume_equation(region, forest, district, species_code, eq_type)

        # --- 3. Prepare ctypes variables for the DLL call ---
        c_regn = c_int(region)
        c_forst = create_string_buffer(f'{forest:<2}'.encode('ascii'))
        c_dist_str = f'{int(district):02d}'
        c_idist = c_int(int(district))
        c_voleq_in = create_string_buffer(f'{voleq_str:<10}'.encode('ascii'))
        c_fiaspcd = c_int(species_code)
        c_dbhob = c_float(dbh)
        c_httot = c_float(height)
        c_mtopp = c_float(merch_top_dib_primary)
        c_mtops = c_float(merch_top_dib_secondary)
        c_stump = c_float(stump_height)

        c_cr = c_float(crown_ratio)
        c_cull = c_float(cull)
        c_decaycd = c_int(decay_cd)
        c_htbrk = c_float(broken_ht)
        c_live = create_string_buffer(f'{live_dead.upper():<1}'.encode('ascii'))

        ctype_char = 'F'
        if calc_type.upper() == 'FIA':
            ctype_char = 'I'
        c_mctype = create_string_buffer(f'{ctype_char:<1}'.encode('ascii'))

        # --- 4. Initialize Output Arrays and other required (but often unused) args ---
        REAL = c_float
        VOL = (REAL * 15)()
        DRYBIO = (REAL * 15)()
        GRNBIO = (REAL * 15)()
        LOGVOL = (REAL * 7 * 20)()
        LOGDIA = (REAL * 21 * 3)()
        LOGLEN = (REAL * 20)()
        BOLHT = (REAL * 21)()
        TLOGS = c_int(0)

        # Set flags for which volumes to calculate
        flags = kwargs.get('flags', {'cutflg': 1, 'cupflg': 1, 'bfpflg': 1})
        c_cutflg = c_int(flags.get('cutflg', 0)) # Total Cubic Vol
        c_cupflg = c_int(flags.get('cupflg', 0)) # Merch Cubic Primary
        c_bfpflg = c_int(flags.get('bfpflg', 0)) # Board Foot Primary
        c_cdpflg = c_int(flags.get('cdpflg', 0)) # Cords Primary
        c_spflg = c_int(flags.get('spflg', 0))  # Secondary Product

        # *** CORRECTED THIS SECTION ***
        # Placeholder variables (floats vs. ints) are now correctly typed
        c_drcob, c_ht1prd, c_ht2prd, c_upsht1, c_upsht2, c_upsd1, c_upsd2 = (c_float(0),) * 7
        c_avgz1, c_avgz2, c_dbtbh, c_btr = (c_float(0),) * 4
        c_nologp, c_nologs = c_float(0), c_float(0)
        c_htbrkd = c_float(0)
        
        # Integer placeholders
        c_htlog, c_htref = c_int(0), c_int(0)
        c_fclass = c_int(0)
        c_httfll, c_ba, c_si = c_int(0), c_int(0), c_int(0)
        # *** END OF CORRECTION ***

        c_conspec = create_string_buffer(b'    ')
        c_prod_vol = create_string_buffer(b'01')
        c_errflag_vol = c_int(0)
        c_httype = create_string_buffer(b'F')

        # --- 5. Call the VOLUMELIBRARY2 function ---
        self.vollib.VOLUMELIBRARY2(
            byref(c_regn), c_forst, len(c_forst)-1, c_voleq_in, len(c_voleq_in)-1,
            byref(c_mtopp), byref(c_mtops), byref(c_stump), byref(c_dbhob),
            byref(c_drcob), c_httype, len(c_httype)-1, byref(c_httot), byref(c_htlog),
            byref(c_ht1prd), byref(c_ht2prd), byref(c_upsht1), byref(c_upsht2),
            byref(c_upsd1), byref(c_upsd2), byref(c_htref), byref(c_avgz1),
            byref(c_avgz2), byref(c_fclass), byref(c_dbtbh), byref(c_btr),
            byref(VOL), byref(LOGVOL), byref(LOGDIA), byref(LOGLEN), byref(BOLHT),
            byref(TLOGS), byref(c_nologp), byref(c_nologs), byref(c_cutflg),
            byref(c_bfpflg), byref(c_cupflg), byref(c_cdpflg), byref(c_spflg),
            c_conspec, len(c_conspec)-1, c_prod_vol, len(c_prod_vol)-1,
            byref(c_httfll), c_live, len(c_live)-1, byref(c_ba), byref(c_si),
            c_mctype, len(c_mctype)-1, byref(c_errflag_vol), byref(c_idist),
            byref(c_htbrk), byref(c_htbrkd), byref(c_fiaspcd), byref(DRYBIO),
            byref(GRNBIO), byref(c_cr), byref(c_cull), byref(c_decaycd)
        )

        # --- 6. Handle errors and return results ---
        if c_errflag_vol.value != 0:
            error_message = VOLLIB_ERROR_CODES.get(c_errflag_vol.value, f"Unknown error code {c_errflag_vol.value} from VOLUMELIBRARY2.")
            raise ValueError(f"{error_message} (Voleq: {voleq_str})")

        return {
            "voleq": voleq_str,
            "VOL": VOL,
            "DRYBIO": DRYBIO,
            "GRNBIO": GRNBIO,
            "BOLHT": BOLHT,
            "LOGDIA": LOGDIA,
            "TLOGS": TLOGS.value
        }

    # --- Public Methods (mimicking VBA functions) ---

    def get_full_report(self, **kwargs):
        """
        Calculates all primary outputs and returns them in a dictionary,
        including volume, biomass, and taper data.
        """
        results = self._call_volume_library(**kwargs)
        
        # Extract taper data into a DataFrame
        heights = [h for h in results["BOLHT"] if h > 0]
        dob = [results["LOGDIA"][i][0] for i, h in enumerate(results["BOLHT"]) if h > 0]
        dib = [results["LOGDIA"][i][1] for i, h in enumerate(results["BOLHT"]) if h > 0]
        taper_df = pd.DataFrame({"height_ft": heights, "dob_in": dob, "dib_in": dib})

        return {
            "voleq": results["voleq"],
            "total_cuft": results["VOL"][0],
            "merch_cuft_primary": results["VOL"][3],
            "merch_board_feet_primary": results["VOL"][1], # Scribner
            "merch_cords_primary": results["VOL"][5],
            "agb_dry_lbs": results["DRYBIO"][0],
            "agb_green_lbs": results["GRNBIO"][0],
            "taper_data": taper_df,
            "number_of_logs": results["TLOGS"]
        }

    def calc_total_cubic_ft(self, **kwargs) -> float:
        """Calculates Total Cubic Foot volume (ground to tip). Corresponds to VOL(1)."""
        kwargs['flags'] = {'cutflg': 1}
        results = self._call_volume_library(**kwargs)
        return results["VOL"][0]

    def calc_merch_cubic_ft(self, **kwargs) -> float:
        """Calculates Gross Merchantable Cubic Foot volume. Corresponds to VOL(4)."""
        kwargs['flags'] = {'cupflg': 1}
        results = self._call_volume_library(**kwargs)
        return results["VOL"][3]

    def calc_merch_board_feet(self, **kwargs) -> float:
        """Calculates Gross Scribner Board Foot volume. Corresponds to VOL(2)."""
        kwargs['flags'] = {'bfpflg': 1}
        results = self._call_volume_library(**kwargs)
        return results["VOL"][1]

    def calc_merch_cords(self, **kwargs) -> float:
        """Calculates Merchantable Cordwood volume. Corresponds to VOL(6)."""
        kwargs['flags'] = {'cdpflg': 1, 'cupflg': 1} # cupflg is often required
        results = self._call_volume_library(**kwargs)
        return results["VOL"][5]

    def _calc_biomass_component(self, dry_or_green: str, index: int, **kwargs) -> float:
        """Helper for all biomass calculations."""
        kwargs['eq_type'] = 'biomass'
        kwargs['flags'] = {'cutflg': 1} # Enable volume calc for biomass dependencies
        results = self._call_volume_library(**kwargs)
        
        if dry_or_green.lower().startswith('d'):
            return results["DRYBIO"][index]
        elif dry_or_green.lower().startswith('g'):
            return results["GRNBIO"][index]
        else:
            raise ValueError("dry_or_green must be 'dry' or 'green'")

    def calc_wt_above_ground_biomass(self, dry_or_green: str, **kwargs) -> float:
        """Calculates above ground biomass (no foliage). Corresponds to DRYBIO/GRNBIO[1]."""
        return self._calc_biomass_component(dry_or_green, 0, **kwargs)

    def calc_wt_total_stem(self, dry_or_green: str, **kwargs) -> float:
        """Calculates total stem biomass (wood + bark). Corresponds to DRYBIO/GRNBIO[2] + [3]."""
        wood = self._calc_biomass_component(dry_or_green, 1, **kwargs)
        bark = self._calc_biomass_component(dry_or_green, 2, **kwargs)
        return wood + bark

    def calc_wt_merch_stem(self, dry_or_green: str, **kwargs) -> float:
        """Calculates merchantable stem biomass (wood + bark). Corresponds to DRYBIO/GRNBIO[6] + [7]."""
        wood = self._calc_biomass_component(dry_or_green, 5, **kwargs)
        bark = self._calc_biomass_component(dry_or_green, 6, **kwargs)
        return wood + bark
        
    def calc_wt_branches(self, dry_or_green: str, **kwargs) -> float:
        """Calculates branch biomass. Corresponds to DRYBIO/GRNBIO[12]."""
        return self._calc_biomass_component(dry_or_green, 11, **kwargs)

    def calc_wt_foliage(self, dry_or_green: str, **kwargs) -> float:
        """Calculates foliage biomass. Corresponds to DRYBIO/GRNBIO[13]."""
        return self._calc_biomass_component(dry_or_green, 12, **kwargs)

    def calc_carbon_content(self, **kwargs) -> float:
        """Calculates carbon content in pounds. Corresponds to DRYBIO[15]."""
        # Carbon is always based on dry weight.
        return self._calc_biomass_component('dry', 14, **kwargs)


# --- EXAMPLE USAGE ---
if __name__ == '__main__':
    # IMPORTANT: Update this path to the location of your vollib DLL file
    # Ensure you use the correct version (vollib.dll for 32-bit Python, vollib64.dll for 64-bit)
    dll_file_path = r"D:/downloaded/forest_analytics/vol-lib-dll-20250701/VolLibDll20250701/vollib64/vollib.dll"

    try:
        # 1. Create an instance of the wrapper
        vol_estimator = VollibWrapper(dll_path=dll_file_path)

        # 2. Define tree parameters in a dictionary
        # Example: Douglas-fir in Region 1, Lolo National Forest (16), Missoula District (01)
        tree_params = {
            "region": 1,
            "forest": "16",
            "district": "01",
            "species_code": 202, # Douglas-fir
            "dbh": 15.8,          # Diameter at Breast Height (inches)
            "height": 88.0,       # Total tree height (feet)
        }
        
        # --- 3. Call the specific methods you need ---
        print("\n--- Individual Calculations ---")
        
        # Calculate total cubic foot volume
        total_cuft = vol_estimator.calc_total_cubic_ft(**tree_params)
        print(f"Total Cubic Foot Volume: {total_cuft:.2f} cuft")

        # Calculate merchantable board feet to a 6" top
        merch_bf = vol_estimator.calc_merch_board_feet(
            **tree_params,
            merch_top_dib_primary=0.0
        )
        print(f"Merchantable Board Feet (to 6\" top): {merch_bf:.0f} BF")

        # Calc total merchantable cubic feet
        merch_cuft = vol_estimator.calc_merch_cubic_ft(
            **tree_params,
            merch_top_dib_primary=0.0
        )
        print(f"Merchantable Cubic Feet (to 6\" top): {merch_cuft:.2f} cuft")

        # Calculate green above-ground biomass
        # Note: Biomass calculations require more FVS-style inputs
        biomass_params = {**tree_params, 'crown_ratio': 40, 'live_dead': 'L'}
        agb_green = vol_estimator.calc_wt_above_ground_biomass(
            dry_or_green='green',
            **biomass_params
        )
        print(f"Above Ground Biomass (Green): {agb_green:.1f} lbs")

        # Calculate carbon content
        carbon = vol_estimator.calc_carbon_content(**biomass_params)
        print(f"Carbon Content: {carbon:.1f} lbs")

        # --- 4. Or, get a full report with taper data ---
        print("\n--- Full Report ---")
        full_report = vol_estimator.get_full_report(
             **tree_params,
             merch_top_dib_primary=6.0
        )
        print(f"Volume Equation Used: {full_report['voleq']}")
        print(f"Total cuft: {full_report['total_cuft']:.2f}")
        print(f"Merch BF: {full_report['merch_board_feet_primary']:.0f}")
        print("\nTaper Profile:")
        print(full_report['taper_data'].to_string(index=False))

    except (OSError, ValueError) as e:
        print(f"\nAn error occurred: {e}")