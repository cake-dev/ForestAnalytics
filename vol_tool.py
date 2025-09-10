import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ctypes import *
import argparse
import sys

# --- 1. Error Codes & Core Wrapper Class ---

VOLLIB_ERROR_CODES = {1:"Invalid Voleq", 2:"Missing Form Class", 3:"DBH < 1.0", 4:"Height < 4.5", 5:"D2H out of range", 6:"Invalid Species Code", 7:"Illegal Primary Log Height", 8:"Illegal Secondary Log Height", 9:"Upper Stem Diameter required", 10:"Illegal Upper Stem Height", 11:"Profile Fit Error", 12:"Tree has > 20 logs", 13:"Top DIB > DBHIB", 14:"Bark Equation Error", 15:"Invalid Biomass Equation", 16:"HT1PRD required for biomass", 17:"HT2PRD required for biomass", 18:"CV4 needed for biomass", 19:"Not a profile model", 20:"Invalid Region", 21:"Invalid Forest", 22:"Invalid District", 23:"Invalid FIA Species Code"}

class VollibWrapper:
    """A comprehensive Python wrapper for the NVEL (VOLLIB) FORTRAN library."""
    def __init__(self, dll_path: str):
        self.dll_path = dll_path
        try:
            self.vollib = windll.LoadLibrary(self.dll_path)
            print(f"✅ Successfully loaded vollib.dll from: {self.dll_path}")
        except OSError as e:
            print(f"❌ Error loading DLL: {e}", file=sys.stderr)
            raise
        self._define_function_prototypes()

    def _define_function_prototypes(self):
        """Sets the argument types (argtypes) for all required DLL functions."""
        # Main calculation engine
        self.vollib.VOLUMELIBRARY2.argtypes=[POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),c_char_p,c_int,POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER((c_float*15)),POINTER((c_float*7*20)),POINTER((c_float*21*3)),POINTER((c_float*20)),POINTER((c_float*21)),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_int),c_char_p,c_int,POINTER(c_int),POINTER(c_int),c_char_p,c_int,POINTER(c_int),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER((c_float*15)),POINTER((c_float*15)),POINTER(c_float),POINTER(c_float),POINTER(c_int)]
        # Equation Lookups
        self.vollib.GETVOLEQ.argtypes=[POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_int)]
        self.vollib.GETNVBEQ2.argtypes=[POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_int),c_char_p,c_int,POINTER(c_int)]
        # Taper/Height functions
        self.vollib.CALCDIA.argtypes=[POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_int)]
        self.vollib.CALCHT2TOPD.argtypes=[POINTER(c_int),c_char_p,c_int,c_char_p,c_int,POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int),POINTER(c_float),POINTER(c_float),POINTER(c_int)]
        # Utility functions
        self.vollib.VERNUM.argtypes=[POINTER(c_int)]
        self.vollib.GETWTFAC.argtypes=[POINTER(c_int),c_char_p,c_int,POINTER(c_int),POINTER(c_float)]

    def _call_volume_library(self, **kwargs):
        """Internal engine that calls VOLUMELIBRARY2 and returns all result arrays."""
        # This function prepares all inputs and calls the main DLL subroutine
        region=kwargs.get('region',1);forest=kwargs.get('forest','01');district=kwargs.get('district','01');species_code=kwargs.get('species_code');dbh=kwargs.get('dbh');height=kwargs.get('height');voleq_override=kwargs.get('voleq_override');merch_top_dib_primary=kwargs.get('merch_top_dib_primary',4.0);merch_top_dib_secondary=kwargs.get('merch_top_dib_secondary',0.0);stump_height=kwargs.get('stump_height',1.0);crown_ratio=kwargs.get('crown_ratio',0.0);cull=kwargs.get('cull',0.0);decay_cd=kwargs.get('decay_cd',0);broken_ht=kwargs.get('broken_ht',0.0);live_dead=kwargs.get('live_dead','L');calc_type=kwargs.get('calc_type','FVS');eq_type=kwargs.get('eq_type','volume')
        voleq_str=voleq_override if voleq_override else self.get_volume_equation(region,forest,district,species_code,eq_type)
        c_regn=c_int(region);c_forst=create_string_buffer(f'{forest:<2}'.encode('ascii'));c_idist=c_int(int(district));c_voleq_in=create_string_buffer(f'{voleq_str:<10}'.encode('ascii'));c_fiaspcd=c_int(species_code);c_dbhob=c_float(dbh);c_httot=c_float(height);c_mtopp=c_float(merch_top_dib_primary);c_mtops=c_float(merch_top_dib_secondary);c_stump=c_float(stump_height);c_cr=c_float(crown_ratio);c_cull=c_float(cull);c_decaycd=c_int(decay_cd);c_htbrk=c_float(broken_ht);c_live=create_string_buffer(f'{live_dead.upper():<1}'.encode('ascii'));c_mctype=create_string_buffer(b'F' if calc_type.upper()!='FIA' else b'I')
        VOL=(c_float*15)();DRYBIO=(c_float*15)();GRNBIO=(c_float*15)();LOGVOL=(c_float*7*20)();LOGDIA=(c_float*21*3)();LOGLEN=(c_float*20)();BOLHT=(c_float*21)();TLOGS=c_int(0)
        flags=kwargs.get('flags',{'cutflg':1,'cupflg':1,'bfpflg':1,'cdpflg':1,'spflg':1})
        c_cutflg=c_int(flags.get('cutflg',0));c_cupflg=c_int(flags.get('cupflg',0));c_bfpflg=c_int(flags.get('bfpflg',0));c_cdpflg=c_int(flags.get('cdpflg',0));c_spflg=c_int(flags.get('spflg',0))
        c_drcob,c_ht1prd,c_ht2prd,c_upsht1,c_upsht2,c_upsd1,c_upsd2=(c_float(0),)*7;c_avgz1,c_avgz2,c_dbtbh,c_btr,c_nologp,c_nologs,c_htbrkd=(c_float(0),)*7;c_htlog,c_htref,c_fclass,c_httfll,c_ba,c_si=(c_int(0),)*6;c_conspec=create_string_buffer(b'    ');c_prod_vol=create_string_buffer(b'01');c_errflag_vol=c_int(0);c_httype=create_string_buffer(b'F')
        self.vollib.VOLUMELIBRARY2(byref(c_regn),c_forst,len(c_forst)-1,c_voleq_in,len(c_voleq_in)-1,byref(c_mtopp),byref(c_mtops),byref(c_stump),byref(c_dbhob),byref(c_drcob),c_httype,len(c_httype)-1,byref(c_httot),byref(c_htlog),byref(c_ht1prd),byref(c_ht2prd),byref(c_upsht1),byref(c_upsht2),byref(c_upsd1),byref(c_upsd2),byref(c_htref),byref(c_avgz1),byref(c_avgz2),byref(c_fclass),byref(c_dbtbh),byref(c_btr),byref(VOL),byref(LOGVOL),byref(LOGDIA),byref(LOGLEN),byref(BOLHT),byref(TLOGS),byref(c_nologp),byref(c_nologs),byref(c_cutflg),byref(c_bfpflg),byref(c_cupflg),byref(c_cdpflg),byref(c_spflg),c_conspec,len(c_conspec)-1,c_prod_vol,len(c_prod_vol)-1,byref(c_httfll),c_live,len(c_live)-1,byref(c_ba),byref(c_si),c_mctype,len(c_mctype)-1,byref(c_errflag_vol),byref(c_idist),byref(c_htbrk),byref(c_htbrkd),byref(c_fiaspcd),byref(DRYBIO),byref(GRNBIO),byref(c_cr),byref(c_cull),byref(c_decaycd))
        if c_errflag_vol.value!=0: raise ValueError(f"{VOLLIB_ERROR_CODES.get(c_errflag_vol.value, f'Unknown error code {c_errflag_vol.value}')} (Voleq: {voleq_str})")
        return {"voleq":voleq_str,"VOL":VOL,"DRYBIO":DRYBIO,"GRNBIO":GRNBIO,"BOLHT":BOLHT,"LOGDIA":LOGDIA,"LOGVOL":LOGVOL,"LOGLEN":LOGLEN,"TLOGS":TLOGS.value}

    # --- Utility Methods ---
    def get_version_number(self) -> int:
        """Returns the version number of the DLL."""
        version = c_int(0)
        self.vollib.VERNUM(byref(version))
        return version.value
    def get_volume_equation(self, region, forest, district, species_code, eq_type='volume') -> str:
        """Looks up the default volume or biomass equation for a given location and species."""
        c_regn=c_int(region);c_forst=create_string_buffer(f'{forest:<2}'.encode('ascii'));c_dist=create_string_buffer(f'{district:<2}'.encode('ascii'));c_spec=c_int(species_code);c_voleq=create_string_buffer(11);c_errflag=c_int(0)
        if eq_type=='volume': self.vollib.GETVOLEQ(byref(c_regn),c_forst,2,c_dist,2,byref(c_spec),create_string_buffer(b'01'),2,c_voleq,10,byref(c_errflag))
        elif eq_type=='biomass': self.vollib.GETNVBEQ2(byref(c_regn),c_forst,2,c_dist,2,byref(c_spec),c_voleq,10,byref(c_errflag))
        if c_errflag.value!=0: raise ValueError(VOLLIB_ERROR_CODES.get(c_errflag.value, f"Eq Lookup Error {c_errflag.value}"))
        return c_voleq.value.decode('ascii').strip()
    def get_weight_factor(self, region, forest, species_code) -> float:
        """Retrieves moisture content conversion factor for a species."""
        c_regn=c_int(region); c_forst=create_string_buffer(f'{forest:<2}'.encode('ascii')); c_spec=c_int(species_code); c_wf=c_float(0)
        self.vollib.GETWTFAC(byref(c_regn), c_forst, 2, byref(c_spec), byref(c_wf))
        return c_wf.value

    # --- Taper, Height, and Diameter Methods ---
    def get_taper_profile_data(self, **kwargs) -> pd.DataFrame:
        """Generates a DataFrame with the tree's taper profile."""
        results = self._call_volume_library(**kwargs)
        heights = [h for h in results["BOLHT"] if h > 0]; dob = [results["LOGDIA"][i][2] for i,h in enumerate(results["BOLHT"]) if h>0]; dib = [results["LOGDIA"][i][1] for i,h in enumerate(results["BOLHT"]) if h>0]
        return pd.DataFrame({"height_ft": heights, "dob_in": dob, "dib_in": dib})
    def plot_taper_profile(self, taper_df: pd.DataFrame, title: str, output_file=None):
        """Uses matplotlib to plot the taper profile from a DataFrame."""
        plt.style.use('seaborn-v0_8-whitegrid');fig,ax=plt.subplots(figsize=(8,10));ax.plot(taper_df['dob_in'],taper_df['height_ft'],label='Outside Bark (DOB)',color='saddlebrown',linewidth=2);ax.plot(taper_df['dib_in'],taper_df['height_ft'],label='Inside Bark (DIB)',color='darkorange',linestyle='--',linewidth=2);ax.set_xlabel("Diameter (inches)",fontsize=12);ax.set_ylabel("Height Above Ground (ft)",fontsize=12);ax.set_title(title,fontsize=14,weight='bold');ax.legend();ax.set_ylim(bottom=0);ax.set_xlim(left=0);fig.tight_layout()
        if output_file:
            try: plt.savefig(output_file); print(f"📈 Plot saved to {output_file}")
            except Exception as e: print(f"❌ Could not save plot: {e}", file=sys.stderr)
        plt.show()
    def calc_dib_at_height(self, height_at, **kwargs) -> float:
        """Calculates diameter inside bark at a specific height on the stem."""
        c_regn=c_int(kwargs['region']);c_forst=create_string_buffer(f"{kwargs['forest']:<2}".encode('ascii'));c_voleq=create_string_buffer(f"{self.get_volume_equation(**kwargs):<10}".encode('ascii'));c_dbhob=c_float(kwargs['dbh']);c_httot=c_float(kwargs['height']);c_htup=c_float(height_at);c_dib=c_float(0);c_dob=c_float(0);c_err=c_int(0)
        c_stump,c_drcob,c_upsht1,c_upsht2,c_upsd1,c_upsd2,c_avgz1,c_avgz2,c_dbtbh,c_btr=(c_float(0),)*10;c_htref,c_fclass=(c_int(0),)*2
        self.vollib.CALCDIA(byref(c_regn),c_forst,2,c_voleq,10,byref(c_stump),byref(c_dbhob),byref(c_drcob),byref(c_httot),byref(c_upsht1),byref(c_upsht2),byref(c_upsd1),byref(c_upsd2),byref(c_htref),byref(c_avgz1),byref(c_avgz2),byref(c_fclass),byref(c_dbtbh),byref(c_btr),byref(c_htup),byref(c_dib),byref(c_dob),byref(c_err))
        if c_err.value!=0: raise ValueError(VOLLIB_ERROR_CODES.get(c_err.value, f"CALCDIA Error {c_err.value}"))
        return c_dib.value
    def calc_height_to_top_dib(self, top_dib, **kwargs) -> float:
        """Calculates the height on the stem to a specific top diameter inside bark."""
        c_regn=c_int(kwargs['region']);c_forst=create_string_buffer(f"{kwargs['forest']:<2}".encode('ascii'));c_voleq=create_string_buffer(f"{self.get_volume_equation(**kwargs):<10}".encode('ascii'));c_dbhob=c_float(kwargs['dbh']);c_httot=c_float(kwargs['height']);c_stemdib=c_float(top_dib);c_stemht=c_float(0);c_err=c_int(0)
        c_ht1prd,c_ht2prd,c_upsht1,c_upsht2,c_upsd1,c_upsd2,c_avgz1,c_avgz2,c_dbtbh,c_btr=(c_float(0),)*10;c_htref,c_fclass=(c_int(0),)*2
        self.vollib.CALCHT2TOPD(byref(c_regn),c_forst,2,c_voleq,10,byref(c_dbhob),byref(c_httot),byref(c_ht1prd),byref(c_ht2prd),byref(c_upsht1),byref(c_upsht2),byref(c_upsd1),byref(c_upsd2),byref(c_avgz1),byref(c_avgz2),byref(c_htref),byref(c_dbtbh),byref(c_btr),byref(c_fclass),byref(c_stemdib),byref(c_stemht),byref(c_err))
        if c_err.value!=0: raise ValueError(VOLLIB_ERROR_CODES.get(c_err.value, f"CALCHT2TOPD Error {c_err.value}"))
        return c_stemht.value

    # --- Volume Calculation Methods ---
    def _get_vol(self, index, **kwargs) -> float: return self._call_volume_library(**kwargs)["VOL"][index]
    def calc_total_cubic_ft(self, **kwargs) -> float: return self._get_vol(0, **kwargs)
    def calc_merch_board_feet(self, **kwargs) -> float: return self._get_vol(1, **kwargs)
    def calc_merch_cubic_ft(self, **kwargs) -> float:
        """Calculates Gross Merchantable Cubic Foot volume. Corresponds to VOL(4)."""
        kwargs['flags'] = {'cupflg': 1}
        results = self._call_volume_library(**kwargs)
        return results["VOL"][3]
    def calc_merch_cords(self, **kwargs) -> float: return self._get_vol(5, **kwargs)
    def calc_cubic_topwood(self, **kwargs) -> float: return self._get_vol(6, **kwargs)
    def calc_international_board_feet(self, **kwargs) -> float: return self._get_vol(9, **kwargs)
    def calc_board_feet_topwood(self, **kwargs) -> float: return self._get_vol(11, **kwargs)
    def calc_cubic_stump(self, **kwargs) -> float: return self._get_vol(13, **kwargs)
    def calc_cubic_tip(self, **kwargs) -> float: return self._get_vol(14, **kwargs)
    
    # --- Biomass Calculation Methods ---
    def _get_bio(self, dry_or_green, index, **kwargs) -> float:
        kwargs['eq_type']='biomass'; results=self._call_volume_library(**kwargs)
        if dry_or_green.lower().startswith('d'): return results["DRYBIO"][index]
        return results["GRNBIO"][index]
    def calc_wt_above_ground_biomass(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 0, **kwargs)
    def calc_wt_total_stem_wood(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 1, **kwargs)
    def calc_wt_total_stem_bark(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 2, **kwargs)
    def calc_wt_stump_wood(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 3, **kwargs)
    def calc_wt_stump_bark(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 4, **kwargs)
    def calc_wt_merch_stem_wood(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 5, **kwargs)
    def calc_wt_merch_stem_bark(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 6, **kwargs)
    def calc_wt_topwood_wood(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 7, **kwargs)
    def calc_wt_topwood_bark(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 8, **kwargs)
    def calc_wt_tip_wood(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 9, **kwargs)
    def calc_wt_tip_bark(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 10, **kwargs)
    def calc_wt_branches(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 11, **kwargs)
    def calc_wt_foliage(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 12, **kwargs)
    def calc_wt_top_and_limbs(self, dry_or_green, **kwargs) -> float: return self._get_bio(dry_or_green, 13, **kwargs)
    def calc_carbon_content(self, **kwargs) -> float: return self._get_bio('dry', 14, **kwargs)

    # --- Log Calculation Methods ---
    def _get_log_val(self, log_num, array_name, idx1, idx2=None):
        if log_num < 1: raise ValueError("log_num must be 1 or greater.")
        results = self._call_volume_library(flags={'bfpflg':1, 'cupflg':1})
        if log_num > results['TLOGS']: raise ValueError(f"Invalid log_num {log_num}. Tree only has {results['TLOGS']} logs.")
        array = results[array_name]
        return array[log_num-1][idx1] if idx2 is None else array[log_num][idx1] # LOGDIA is 1-based on the first index in VBA
    def calc_log_volume_cubic(self, log_num, **kwargs) -> float: return self._get_log_val(log_num, "LOGVOL", 3, **kwargs)
    def calc_log_volume_board_feet(self, log_num, **kwargs) -> float: return self._get_log_val(log_num, "LOGVOL", 0, **kwargs)
    def calc_log_length(self, log_num, **kwargs) -> float: return self._get_log_val(log_num, "LOGLEN", None, **kwargs)
    def calc_log_dib(self, log_num, **kwargs) -> float: return self._get_log_val(log_num, "LOGDIA", 1, **kwargs)
    def calc_log_scaling_dib(self, log_num, **kwargs) -> float: return self._get_log_val(log_num, "LOGDIA", 0, **kwargs)

def main():
    """Main function to parse command-line arguments and run calculations."""
    parser = argparse.ArgumentParser(description="A command-line tool for the National Volume Estimator Library (VOLLIB).", formatter_class=argparse.RawTextHelpFormatter)
    parent_parser = argparse.ArgumentParser(add_help=False)
    # parent_parser.add_argument('--dll-path', required=True, help="Path to the vollib64.dll file.")
    parent_parser.add_argument('--species-code', type=int, required=True, help="FIA species code (e.g., 202 for Douglas-fir).")
    parent_parser.add_argument('--dbh', type=float, required=True, help="Diameter at breast height (inches).")
    parent_parser.add_argument('--height', type=float, required=True, help="Total tree height (feet).")
    parent_parser.add_argument('--region', type=int, default=1, help="FS region (default: 1).")
    parent_parser.add_argument('--forest', type=str, default='16', help="FS forest code (default: '16').")
    parent_parser.add_argument('--district', type=str, default='01', help="FS district code (default: '01').")
    subparsers = parser.add_subparsers(dest='command', required=True, title="Available Commands")
    
    # --- Build CLI for all functions ---
    # Utility
    subparsers.add_parser('version', parents=[parent_parser], help="Get the DLL version number.")
    # Volume
    subparsers.add_parser('total_cuft', parents=[parent_parser], help="Calculate total cubic foot volume.")
    p_mbf=subparsers.add_parser('merch_bf', parents=[parent_parser], help="Calculate Scribner board foot volume."); p_mbf.add_argument('--top-dib',type=float,default=6.0,help="Merch top DIB (default: 6.0).")
    p_mcf=subparsers.add_parser('merch_cuft', parents=[parent_parser], help="Calculate merch. cubic foot volume."); p_mcf.add_argument('--top-dib',type=float,default=4.0,help="Merch top DIB (default: 4.0).")
    p_int=subparsers.add_parser('intl_bf', parents=[parent_parser], help="Calculate International 1/4\" board foot volume."); p_int.add_argument('--top-dib',type=float,default=6.0,help="Merch top DIB (default: 6.0).")
    # Biomass
    p_agb=subparsers.add_parser('agb', parents=[parent_parser], help="Calculate above-ground biomass."); p_agb.add_argument('dry_or_green', choices=['dry','green']); p_agb.add_argument('--crown-ratio',type=int,default=40)
    p_c=subparsers.add_parser('carbon', parents=[parent_parser], help="Calculate total carbon content in pounds."); p_c.add_argument('--crown-ratio',type=int,default=40)
    # Taper
    p_dib=subparsers.add_parser('dib_at_ht', parents=[parent_parser], help="Calculate DIB at a specific height."); p_dib.add_argument('height_at',type=float,help="Height on stem (ft).")
    p_ht=subparsers.add_parser('ht_to_dib', parents=[parent_parser], help="Calculate height to a specific DIB."); p_ht.add_argument('top_dib',type=float,help="Top DIB (inches).")
    p_plot=subparsers.add_parser('plot_taper', parents=[parent_parser], help="Generate and display a taper profile plot."); p_plot.add_argument('--output-file', type=str, help="Optional path to save the plot image.")
    # Logs
    p_logc=subparsers.add_parser('log_cuft', parents=[parent_parser], help="Get cubic foot volume for a specific log."); p_logc.add_argument('log_num',type=int)
    p_logb=subparsers.add_parser('log_bf', parents=[parent_parser], help="Get board foot volume for a specific log."); p_logb.add_argument('log_num',type=int)
    p_logl=subparsers.add_parser('log_len', parents=[parent_parser], help="Get the length of a specific log."); p_logl.add_argument('log_num',type=int)
    p_logd=subparsers.add_parser('log_dib', parents=[parent_parser], help="Get the small-end DIB of a specific log."); p_logd.add_argument('log_num',type=int)

    args = parser.parse_args()
    params = vars(args)
    if 'top_dib' in params: params['merch_top_dib_primary'] = params.pop('top_dib')
        
    try:
        vol_estimator = VollibWrapper(dll_path=r"D:/downloaded/forest_analytics/vol-lib-dll-20250701/VolLibDll20250701/vollib64/vollib.dll")
        cmd = args.command
        
        # --- Execute Command ---
        if cmd == 'version': result = f"{vol_estimator.get_version_number()}"
        elif cmd == 'total_cuft': result = f"{vol_estimator.calc_total_cubic_ft(**params):.2f} total cubic feet"
        elif cmd == 'merch_bf': result = f"{vol_estimator.calc_merch_board_feet(**params):.0f} Scribner board feet"
        elif cmd == 'merch_cuft': result = f"{vol_estimator.calc_merch_cubic_ft(**params):.2f} merchantable cubic feet"
        elif cmd == 'intl_bf': result = f"{vol_estimator.calc_international_board_feet(**params):.0f} International 1/4\" board feet"
        elif cmd == 'agb': result = f"{vol_estimator.calc_wt_above_ground_biomass(**params):.1f} lbs ({args.dry_or_green})"
        elif cmd == 'carbon': result = f"{vol_estimator.calc_carbon_content(**params):.1f} lbs of carbon"
        elif cmd == 'dib_at_ht': result = f"{vol_estimator.calc_dib_at_height(**params):.2f} inches DIB at {args.height_at} ft"
        elif cmd == 'ht_to_dib': result = f"{vol_estimator.calc_height_to_top_dib(**params):.1f} ft is height to {args.top_dib}\" DIB"
        elif cmd == 'log_cuft': result = f"Log {args.log_num}: {vol_estimator.calc_log_volume_cubic(**params):.2f} cubic feet"
        elif cmd == 'log_bf': result = f"Log {args.log_num}: {vol_estimator.calc_log_volume_board_feet(**params):.0f} board feet"
        elif cmd == 'log_len': result = f"Log {args.log_num}: {vol_estimator.calc_log_length(**params):.1f} feet long"
        elif cmd == 'log_dib': result = f"Log {args.log_num}: {vol_estimator.calc_log_dib(**params):.1f} inches small-end DIB"
        elif cmd == 'plot_taper':
            taper_data=vol_estimator.get_taper_profile_data(**params);eq_params={'region':args.region,'forest':args.forest,'district':args.district,'species_code':args.species_code};voleq=vol_estimator.get_volume_equation(**eq_params);title=f"Taper Profile (Voleq: {voleq})\nDBH={args.dbh}\", Ht={args.height}\"";vol_estimator.plot_taper_profile(taper_data,title=title,output_file=args.output_file);return
        print(f"\nResult: {result}")

    except (OSError, ValueError, TypeError) as e:
        print(f"\n❌ ERROR: {e}", file=sys.stderr); sys.exit(1)

if __name__ == '__main__':
    main()