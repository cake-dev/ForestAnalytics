# Forest Inventory Analysis Using VOLLIB - Assignment 1
# NVEL Volume and Taper Estimation with Lubrecht Plot 9 Data

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from vol_eqs import VollibWrapper
import warnings
warnings.filterwarnings('ignore')

# Set up plotting style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10

# IMPORTANT: Update this path to your vollib DLL location
DLL_PATH = r"D:/downloaded/forest_analytics/vol-lib-dll-20250701/VolLibDll20250701/vollib64/vollib.dll"

print("=== Forest Inventory Analysis Using VOLLIB ===")
print("Analysis of Lubrecht Experimental Forest Plot 9\n")

# Initialize the VOLLIB wrapper
try:
    vol_estimator = VollibWrapper(dll_path=DLL_PATH)
    print("✓ VOLLIB DLL loaded successfully")
except Exception as e:
    print(f"✗ Error loading VOLLIB DLL: {e}")
    print("Please update the DLL_PATH variable to match your system")
    raise

# Load the plot data from Excel file
try:
    plot_data = pd.read_excel('lef_plot9.xlsx', sheet_name=0)  # First sheet
    species_codes = pd.read_csv('fia_treenames.csv')
    plot_data = plot_data.merge(species_codes[['FIA Code', 'Common Name']], left_on='SPCD', right_on='FIA Code', how='left')
    plot_data.rename(columns={'Common Name': 'SpeciesName'}, inplace=True)
    print("✓ Plot data loaded successfully")
    print(f"Dataset contains {len(plot_data)} trees")
except FileNotFoundError:
    print("✗ Excel file 'lef_plot9.xlsx' not found")
    print("Creating sample data for demonstration...")
    
    # Sample data structure
    plot_data = pd.DataFrame({
        'TreeID': range(1, 31), 'PlotID': [9] * 30,
        'SPCD': [202, 202, 93, 122, 93, 202, 42] * 4 + [202, 93],
        'SpeciesName': ['Douglas-fir', 'Douglas-fir', 'Engelmann spruce', 'Lodgepole pine', 'Engelmann spruce', 'Douglas-fir', 'Subalpine fir'] * 4 + ['Douglas-fir', 'Engelmann spruce'],
        'DBH': np.random.normal(12, 4, 30), 'TotHeight': np.random.normal(60, 15, 30),
        'Status': ['Live'] * 28 + ['Dead'] * 2,
    })
    plot_data['DBH'] = np.maximum(plot_data['DBH'], 1.0)
    plot_data['TotHeight'] = np.maximum(plot_data['TotHeight'], 10.0)

# ===================================================================
# TASK 1: Analysis of Tree 24 - Volume Calculations
# ===================================================================

print("\n" + "="*60)
print("TASK 1: Analysis of Tree 24")
print("="*60)

if 24 in plot_data['TreeID'].values:
    tree24 = plot_data[plot_data['TreeID'] == 24].iloc[0]
else:
    tree24 = plot_data[plot_data['SPCD'] == 202].iloc[0]
    print(f"Using TreeID {tree24['TreeID']} as a stand-in for Tree 24 (Douglas-fir)")

print(f"\nTree 24 Characteristics:")
print(f"Species: {tree24['SpeciesName']} (Code: {tree24['SPCD']})")
print(f"DBH: {tree24['DBH']:.1f} inches")
print(f"Total Height: {tree24['TotHeight']:.1f} feet")

tree_params = {
    "region": 1, "forest": "16", "district": "01",
    "species_code": int(tree24['SPCD']),
    "dbh": tree24['DBH'], "height": tree24['TotHeight'],
}

try:
    total_cuft = vol_estimator.calc_total_cubic_ft(**tree_params)
    merch_cuft = vol_estimator.calc_merch_cubic_ft(**tree_params, merch_top_dib_primary=4.0)
    
    basal_area = 0.005454154 * (tree24['DBH'] ** 2)
    
    print(f"\nVolume Calculations for Tree 24:")
    print(f"Basal Area: {basal_area:.3f} sq ft")
    print(f"Total Cubic Volume: {total_cuft:.2f} cu ft")
    print(f"Merchantable Volume (4\" top): {merch_cuft:.2f} cu ft")
    
    tree_params_temp = tree_params.copy()
    # remove dbh and height for equation retrieval
    tree_params_temp.pop('dbh')
    tree_params_temp.pop('height')
    voleq = vol_estimator._get_volume_equation(**tree_params_temp)
    print(f"\nVolume Equation Used: {voleq}")
    
except Exception as e:
    print(f"Error in volume calculations: {e}")
    total_cuft = merch_cuft = 0
    basal_area = 0.005454154 * (tree24['DBH'] ** 2)

print(f"\nEquation Description:")
print("The VOLLIB equations are from the National Volume Estimator Library (NVEL),")
print("which uses region, forest, and species-specific equations developed from")
print("extensive field measurements to predict tree volumes from DBH and height.")

# ===================================================================
# TASK 2: Taper Profile and Girard Form Class
# ===================================================================

print("\n" + "="*60)
print("TASK 2: Taper Profile and Girard Form Class")
print("="*60)

try:
    # *** ADJUSTED: Call get_full_report and extract the taper_data DataFrame ***
    full_report = vol_estimator.get_full_report(**tree_params)
    taper_data = full_report['taper_data']
    
    if not taper_data.empty:
        print(f"\nTaper Profile Data (first 10 points):")
        print(taper_data.head(10).to_string(index=False))
        
        plt.figure(figsize=(10, 6))
        plt.plot(taper_data['height_ft'], taper_data['dob_in'], 
                 'b-', linewidth=2, label='Outside Bark (DOB)')
        plt.plot(taper_data['height_ft'], taper_data['dib_in'], 
                 'r--', linewidth=2, label='Inside Bark (DIB)')
        plt.xlabel('Height Above Ground (ft)')
        plt.ylabel('Diameter (inches)')
        plt.title(f'Taper Profile for Tree 24\nDBH: {tree24["DBH"]:.1f}" Height: {tree24["TotHeight"]:.1f}\'')
        plt.legend()
        plt.show()
        
        # Calculate Girard Form Class: (DIB at 17.3 ft / DBH) * 100
        if 17.3 <= taper_data['height_ft'].max():
            dib_17_3 = np.interp(17.3, taper_data['height_ft'], taper_data['dib_in'])
            gfc = (dib_17_3 / tree24['DBH']) * 100
            print(f"\nGirard Form Class: {gfc:.1f}")
        else:
            print(f"\nTree too short for standard GFC calculation (17.3 ft)")
            
except Exception as e:
    print(f"Error generating taper profile: {e}")

# ===================================================================
# TASK 3: Species Comparison (Douglas-fir vs Engelmann Spruce)
# ===================================================================

print("\n" + "="*60)
print("TASK 3: Species Comparison")
print("="*60)

print("Comparing Tree 24 if it were an Engelmann spruce instead of a Douglas-fir.")

es_params = tree_params.copy()
es_params['species_code'] = 93  # Engelmann Spruce FIA code

try:
    es_total_cuft = vol_estimator.calc_total_cubic_ft(**es_params)
    es_merch_cuft = vol_estimator.calc_merch_cubic_ft(**es_params, merch_top_dib_primary=4.0)
    
    comparison_df = pd.DataFrame({
        'Metric': ['Total Volume (cu ft)', 'Merch Volume (cu ft)'],
        'Douglas-fir': [f"{total_cuft:.2f}", f"{merch_cuft:.2f}"],
        'Engelmann spruce': [f"{es_total_cuft:.2f}", f"{es_merch_cuft:.2f}"]
    })
    print("\nVolume Comparison:")
    print(comparison_df.to_string(index=False))

    # *** ADJUSTED: Use get_full_report for both species' taper profiles ***
    df_taper = vol_estimator.get_full_report(**tree_params)['taper_data']
    es_taper = vol_estimator.get_full_report(**es_params)['taper_data']
    
    if not df_taper.empty and not es_taper.empty:
        plt.figure(figsize=(10, 6))
        plt.plot(df_taper['height_ft'], df_taper['dob_in'], 
                 'g-', linewidth=2, label='Douglas-fir')
        plt.plot(es_taper['height_ft'], es_taper['dob_in'], 
                 'b--', linewidth=2, label='Engelmann spruce')
        plt.xlabel('Height Above Ground (ft)')
        plt.ylabel('Diameter Outside Bark (inches)')
        plt.title('Taper Profile Comparison\n(Same DBH and Height, Different Species)')
        plt.legend()
        plt.show()

    print(f"\nSpecies Differences Explanation:")
    print("Different tree species have different shapes (allometry). Engelmann spruce")
    print("is typically less tapered (more excurrent) than Douglas-fir. The NVEL")
    print("equations capture these species-specific growth patterns, often resulting")
    print("in higher volume estimates for spruce compared to a fir of the same DBH and height.")
    
except Exception as e:
    print(f"Error in species comparison: {e}")

# ===================================================================
# TASK 4: Plot-Level Analysis
# ===================================================================

print("\n" + "="*60)
print("TASK 4: Plot-Level Analysis")
print("="*60)

live_trees = plot_data[(plot_data['Status'] == 1) & (plot_data['DBH'] >= 5.0)].copy()

print(f"Analyzing {len(live_trees)} live trees ≥ 5\" DBH")

def calculate_tree_volumes(row):
    """Calculates volumes and BA for a single tree row."""
    try:
        params = {
            "region": 1, "forest": "16", "district": "01",
            "species_code": int(row['SPCD']),
            "dbh": row['DBH'], "height": row['TotHeight'],
        }
        params_temp = params.copy()
        # remove dbh and height for equation retrieval
        params_temp.pop('dbh', None)
        params_temp.pop('height', None)
        voleq = vol_estimator._get_volume_equation(**params_temp)
        total_vol = vol_estimator.calc_total_cubic_ft(**params)
        merch_vol = vol_estimator.calc_merch_cubic_ft(**params, merch_top_dib_primary=4.0)
        basal_area = 0.005454154 * (row['DBH'] ** 2)
        
        return pd.Series({
            'basalAreasqft': basal_area, 'totVol': total_vol,
            'merchVol': merch_vol, 'volEQnum': voleq
        })
    except Exception:
        basal_area = 0.005454154 * (row['DBH'] ** 2)
        return pd.Series({
            'basalAreasqft': basal_area, 'totVol': 0,
            'merchVol': 0, 'volEQnum': 'Error'
        })

print("Calculating volumes for all trees...")
volume_results = live_trees.apply(calculate_tree_volumes, axis=1)
live_trees = pd.concat([live_trees, volume_results], axis=1)

display_cols = ['TreeID', 'SpeciesName', 'DBH', 'TotHeight', 
                'basalAreasqft', 'totVol', 'merchVol', 'volEQnum']
print(f"\nTree-Level Results:")
print(live_trees[display_cols].round(2).to_string())

# ===================================================================
# TASK 5: Plot Aggregation and Histogram
# ===================================================================

print("\n" + "="*60)
print("TASK 5: Plot Aggregation and Volume Histogram")
print("="*60)

plot_size_acres = 1/9.98  # For trees ≥ 5" DBH
expansion_factor = 1 / plot_size_acres

total_ba_per_acre = live_trees['basalAreasqft'].sum() * expansion_factor
total_vol_per_acre = live_trees['totVol'].sum() * expansion_factor
merch_vol_per_acre = live_trees['merchVol'].sum() * expansion_factor


print(f"Plot Summary (Per-Acre Basis):")
print(f"Plot Size: {plot_size_acres:.4f} acres (expansion factor: {expansion_factor:.1f})")
print(f"Total Basal Area: {total_ba_per_acre} sq ft/acre")
print(f"Total Volume: {total_vol_per_acre:.0f} cu ft/acre")
print(f"Merchantable Volume: {merch_vol_per_acre:.0f} cu ft/acre")

# Create DBH classes for histogram
dbh_min = 5
dbh_max = int(np.ceil(live_trees['DBH'].max()))
bin_width = 2
bins = list(range(dbh_min, dbh_max + bin_width + 1, bin_width))
labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins)-1)]

live_trees['DBH_Class'] = pd.cut(live_trees['DBH'], bins=bins, labels=labels, right=False)
live_trees['merchVol_per_acre'] = live_trees['merchVol'] * expansion_factor

volume_by_class_species = live_trees.pivot_table(
    index='DBH_Class', columns='SpeciesName', 
    values='merchVol_per_acre', aggfunc='sum', fill_value=0
)

print(f"\nMerchantable Volume by DBH Class and Species (cu ft/acre):")
print(volume_by_class_species.round(1))

volume_by_class_species.plot(
    kind='bar', stacked=True, figsize=(12, 8),
    colormap='viridis', edgecolor='black', linewidth=0.5)

plt.title('Merchantable Volume per Acre by DBH Class and Species\nLubrecht Experimental Forest - Plot 9', 
          fontsize=14, pad=20)
plt.xlabel('DBH Class (inches)', fontsize=12)
plt.ylabel('Merchantable Volume (cu ft/acre)', fontsize=12)
plt.xticks(rotation=45)
plt.legend(title='Species')
plt.tight_layout()
plt.show()