################################################################################
#### DELAWARE RIVER BASIN PFAS APP - COMPLETE 1:1 PORT
################################################################################

import streamlit as st
import pandas as pd
import geopandas as gpd
import numpy as np
import folium
from folium import plugins
from streamlit_folium import st_folium
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import warnings
import base64
from io import BytesIO
warnings.filterwarnings('ignore')
from shapely.geometry import Point
import branca.colormap as cm

# Page configuration
st.set_page_config(
    page_title="PFAS in the Delaware River Basin",
    page_icon="🌊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Clean, professional CSS
st.markdown("""
<style>
    .main { background-color: #ffffff; }
    [data-testid="stSidebar"] { background-color: #1a3a52; }
    [data-testid="stSidebar"] .stMarkdown { color: #ffffff; }
    [data-testid="stSidebar"] label { color: #ffffff !important; }
    [data-testid="stSidebar"] .stSelectbox label,
    [data-testid="stSidebar"] .stRadio label,
    [data-testid="stSidebar"] .stSlider label { color: #ffffff !important; }
    .main-header {
        background-color: #1a3a52;
        padding: 2rem 2rem 1.5rem 2rem;
        border-radius: 8px;
        margin-bottom: 2rem;
        border-left: 5px solid #2e7aaf;
    }
    .main-header h1 {
        color: #ffffff;
        font-size: 2.2rem;
        font-weight: 600;
        margin: 0;
        letter-spacing: -0.5px;
    }
    .main-header p {
        color: #b8d4e8;
        font-size: 1rem;
        margin: 0.5rem 0 0 0;
    }
    [data-testid="stMetricValue"] {
        font-size: 1.8rem;
        font-weight: 600;
        color: #1a3a52;
    }
    [data-testid="stMetricLabel"] {
        font-size: 0.9rem;
        color: #666666;
        font-weight: 500;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 4px;
        background-color: #f8f9fa;
        padding: 0.5rem;
        border-radius: 6px;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: transparent;
        border-radius: 4px;
        padding: 0.5rem 1.5rem;
        font-weight: 500;
        color: #495057;
        border: none;
    }
    .stTabs [aria-selected="true"] {
        background-color: #1a3a52;
        color: #ffffff !important;
    }
    .info-card {
        background-color: #f0f7fc;
        border-left: 4px solid #2e7aaf;
        padding: 1.5rem;
        border-radius: 6px;
        margin: 1rem 0;
    }
    .info-card h3 {
        color: #1a3a52;
        margin: 0 0 0.5rem 0;
        font-size: 1.3rem;
        font-weight: 600;
    }
    .info-card p {
        color: #495057;
        margin: 0;
        line-height: 1.6;
    }
    .note-box {
        background-color: #fff8e1;
        border-left: 4px solid #ffc107;
        padding: 1rem;
        border-radius: 4px;
        margin: 1rem 0;
        color: #856404;
    }
    .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        padding: 1rem;
        border-radius: 4px;
        margin: 1rem 0;
        color: #856404;
    }
    .footer {
        background-color: #f8f9fa;
        border-top: 2px solid #1a3a52;
        padding: 2rem;
        border-radius: 6px;
        margin-top: 3rem;
        text-align: center;
        color: #495057;
    }
</style>
""", unsafe_allow_html=True)

# Cache data loading
@st.cache_data
def load_data():
    """Load and preprocess all PFAS data"""

    def simplify_data(df, tissue=False):
        """Simplify and aggregate data"""

        def map_agency(org_name):
            if pd.isna(org_name):
                return "Other"
            org_name = str(org_name)
            if "DRBC" in org_name or "Delaware River Basin Commission" in org_name:
                return "DRBC"
            elif "USGS" in org_name or "U.S. Geological Survey" in org_name:
                return "USGS"
            elif "USEPA" in org_name:
                return "USEPA"
            elif "NJDEP" in org_name or "New Jersey Department of Environmental Protection" in org_name:
                return "NJDEP"
            else:
                return "Other"

        if tissue:
            simple_df = pd.DataFrame({
                'agency': df['OrganizationFormalName'].apply(map_agency),
                'loc': df['Location.Name'],
                'lat': df['Latitude'].round(3),
                'lon': df['Longitude'].round(3),
                'yr': df['Year'],
                'conc': df['Result.Measure.Value..ppt.'],
                'chem': df['PFAS.Chemical.Name'],
                'abbrev': df['PFAS.abbrev'],
                'group': df['PFAS.group'],
                'species': df['Fish.Species.Common']
            })
            cols = ['agency', 'loc', 'yr', 'conc', 'chem', 'abbrev', 'group', 'species']
            simple_df = simple_df.drop_duplicates(subset=cols)
            raw = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'abbrev', 'group', 'agency', 'species'])['conc'].mean().reset_index()
            grouped = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'group', 'agency', 'species'])['conc'].sum().reset_index()
            all_pfas = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'agency', 'species'])['conc'].sum().reset_index()
        else:
            simple_df = pd.DataFrame({
                'agency': df['OrganizationFormalName'].apply(map_agency),
                'loc': df['Location.Name'],
                'lat': df['Latitude'].round(3),
                'lon': df['Longitude'].round(3),
                'yr': df['Year'],
                'conc': df['Result.Measure.Value..ppt.'],
                'chem': df['PFAS.Chemical.Name'],
                'abbrev': df['PFAS.abbrev'],
                'group': df['PFAS.group']
            })
            cols = ['agency', 'loc', 'yr', 'conc', 'chem', 'abbrev', 'group']
            simple_df = simple_df.drop_duplicates(subset=cols)
            raw = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'abbrev', 'group', 'agency'])['conc'].mean().reset_index()
            grouped = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'group', 'agency'])['conc'].sum().reset_index()
            all_pfas = simple_df.groupby(['loc', 'lat', 'lon', 'yr', 'agency'])['conc'].sum().reset_index()

        return {'raw': raw, 'grouped': grouped, 'all': all_pfas}

    def clean_dataframe(df):
        """Remove infinity and invalid values from dataframe"""
        # Replace inf/-inf with NaN
        df = df.replace([np.inf, -np.inf], np.nan)

        # Clean all numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            # Replace any remaining non-finite values with NaN
            df[col] = df[col].where(np.isfinite(df[col]), np.nan)

        # Ensure critical columns are finite
        if 'Latitude' in df.columns:
            df = df[np.isfinite(df['Latitude']) | df['Latitude'].isna()]
        if 'Longitude' in df.columns:
            df = df[np.isfinite(df['Longitude']) | df['Longitude'].isna()]
        if 'Result.Measure.Value..ppt.' in df.columns:
            df = df[np.isfinite(df['Result.Measure.Value..ppt.']) | df['Result.Measure.Value..ppt.'].isna()]

        return df

    data_path = Path(__file__).parent
    water_df = pd.read_csv(data_path / "PFAS_water_data_DRB_20260112.csv")
    gw_df = pd.read_csv(data_path / "PFAS_ground_water_data_DRB_20260112.csv")
    sed_df = pd.read_csv(data_path / "PFAS_sediment_data_DRB_20260112.csv")
    tissue_df = pd.read_csv(data_path / "PFAS_tissue_data_DRB_20260112.csv")

    # Clean all dataframes immediately after loading
    water_df = clean_dataframe(water_df)
    gw_df = clean_dataframe(gw_df)
    sed_df = clean_dataframe(sed_df)
    tissue_df = clean_dataframe(tissue_df)

    sed_df['Result.Measure.Value..ppt.'] = sed_df['Result.Measure.Value..ppt.'] / 1000
    tissue_df['Result.Measure.Value..ppt.'] = tissue_df['Result.Measure.Value..ppt.'] / 1000

    return {
        'Water': simplify_data(water_df),
        'GW': simplify_data(gw_df),
        'Sediment': simplify_data(sed_df),
        'Tissue': simplify_data(tissue_df, tissue=True)
    }

@st.cache_data
def load_spatial_data():
    """Load spatial/shapefile data"""
    data_path = Path(__file__).parent
    drbbnd = gpd.read_file(data_path / "drb_bnd_arc.shp").to_crs(epsg=4326)
    huc = gpd.read_file(data_path / "drbhuc12.shp").to_crs(epsg=4326)
    rm = gpd.read_file(data_path / "delrivRM10th.shp").to_crs(epsg=4326)
    delrivbay = gpd.read_file(data_path / "delrivbay.shp").to_crs(epsg=4326)

    return {'drbbnd': drbbnd, 'huc': huc, 'rm': rm, 'delrivbay': delrivbay}

# Legend configuration
LEGEND_CONFIG = {
    'Water': {
        'breaks': [0, 1, 2, 5, 10, 25, 50, 100, 1000, 999999],
        'labels': ["B.D. - 1", ">1 - 2", ">2 - 5", ">5 - 10", ">10 - 25",
                   ">25 - 50", ">50 - 100", ">100 - 1000", ">1000"],
        'colors': ['#440154', '#443983', '#31688e', '#21918c', '#35b779', '#90d743', '#fde724']
    },
    'GW': {
        'breaks': [0, 1, 2, 5, 10, 25, 50, 100, 1000, 999999],
        'labels': ["B.D. - 1", ">1 - 2", ">2 - 5", ">5 - 10", ">10 - 25",
                   ">25 - 50", ">50 - 100", ">100 - 1000", ">1000"],
        'colors': ['#440154', '#443983', '#31688e', '#21918c', '#35b779', '#90d743', '#fde724']
    },
    'Sediment': {
        'breaks': [0, 0.01, 0.05, 0.10, 0.25, 0.50, 1.00, 2.50, 5.00, 999999],
        'labels': ["B.D. - 0.01", ">0.01 - 0.05", ">0.05 - 0.10", ">0.10 - 0.25",
                   ">0.25 - 0.50", ">0.50 - 1.00", ">1.00 - 2.50", ">2.50 - 5.00", ">5.00"],
        'colors': ['#440154', '#443983', '#31688e', '#21918c', '#35b779', '#90d743', '#fde724']
    },
    'Tissue': {
        'breaks': [0, 0.10, 0.50, 1.00, 5.00, 10.0, 25.0, 50.0, 100.0, 999999],
        'labels': ["B.D. - 0.10", ">0.10 - 0.50", ">0.50 - 1.00", ">1.00 - 5.00",
                   ">5.00 - 10.0", ">10.0 - 25.0", ">25.0 - 50.0", ">50.0 - 100.0", ">100.0"],
        'colors': ['#440154', '#443983', '#31688e', '#21918c', '#35b779', '#90d743', '#fde724']
    }
}

def get_color_for_concentration(conc, dataset_key):
    """Get color based on concentration"""
    breaks = LEGEND_CONFIG[dataset_key]['breaks']
    colors = LEGEND_CONFIG[dataset_key]['colors']

    for i in range(len(breaks) - 1):
        if conc >= breaks[i] and conc < breaks[i + 1]:
            return colors[min(i, len(colors) - 1)]
    return colors[-1]

def get_marker_radius(conc, dataset_key):
    """Get marker radius based on concentration"""
    breaks = LEGEND_CONFIG[dataset_key]['breaks']

    for i in range(len(breaks) - 1):
        if conc >= breaks[i] and conc < breaks[i + 1]:
            # Map to radius 3-12
            return 3 + (i * 1.2)
    return 12

def create_pie_chart_svg(compound_data, group_data):
    """Create SVG pie charts for compound and group breakdown"""

    def make_pie_svg(data, label_col):
        if data.empty:
            return ""

        # Calculate fractions
        total = data['conc'].sum()
        if total == 0:
            return ""

        data = data.copy()
        data['frac'] = data['conc'] / total
        data = data[data['frac'] >= 0.01]  # Only show >= 1%

        # Generate colors
        n = len(data)
        colors = px.colors.sample_colorscale("Portland", [i/(n-1) if n > 1 else 0 for i in range(n)])

        # SVG pie chart
        svg_parts = ['<svg width="140" height="140" viewBox="0 0 140 140">']

        start_angle = 0
        for idx, row in data.iterrows():
            angle = row['frac'] * 360
            end_angle = start_angle + angle

            # Calculate path
            large_arc = 1 if angle > 180 else 0
            start_rad = np.radians(start_angle - 90)
            end_rad = np.radians(end_angle - 90)

            x1 = 70 + 50 * np.cos(start_rad)
            y1 = 70 + 50 * np.sin(start_rad)
            x2 = 70 + 50 * np.cos(end_rad)
            y2 = 70 + 50 * np.sin(end_rad)

            path = f'M 70,70 L {x1},{y1} A 50,50 0 {large_arc},1 {x2},{y2} Z'
            svg_parts.append(f'<path d="{path}" fill="{colors[idx % len(colors)]}" stroke="white" stroke-width="1"/>')

            # Label
            mid_angle = start_angle + angle/2
            mid_rad = np.radians(mid_angle - 90)
            label_x = 70 + 30 * np.cos(mid_rad)
            label_y = 70 + 30 * np.sin(mid_rad)

            label_text = f"{row[label_col]} ({row['frac']*100:.0f}%)"
            svg_parts.append(f'<text x="{label_x}" y="{label_y}" font-size="8" text-anchor="middle" fill="black">{label_text}</text>')

            start_angle = end_angle

        svg_parts.append('</svg>')
        return ''.join(svg_parts)

    compound_svg = make_pie_svg(compound_data, 'abbrev')
    group_svg = make_pie_svg(group_data, 'group')

    if not compound_svg and not group_svg:
        return ""

    return f"""
    <div style="width: 300px;">
        <div><b>By Compound:</b><br>{compound_svg}</div>
        <div style="margin-top:10px;"><b>By Group:</b><br>{group_svg}</div>
    </div>
    """

def create_detailed_popup(lat, lon, full_data, data_type, dataset_key, filter_name):
    """Create detailed popup with all sample information"""

    # Get all data for this location
    loc_data = full_data[(full_data['lat'] == lat) & (full_data['lon'] == lon)].copy()

    if loc_data.empty:
        return "<p>No data</p>"

    unit = "ng/g" if dataset_key in ["Sediment", "Tissue"] else "ng/L"

    # Sort by year
    loc_data = loc_data.sort_values('yr')

    # Header
    html_parts = [f"<div style='font-family: Arial; max-width: 350px; max-height: 400px; overflow-y: auto;'>"]
    html_parts.append(f"<b>{filter_name}</b><br>")
    html_parts.append(f"<b>Location:</b> {loc_data.iloc[0].get('loc', 'Unknown')}<br>")
    html_parts.append("<hr>")

    # If ΣPFAS, show pie charts
    if data_type == "all":
        # Get compound and group breakdown
        compound_data = loc_data.groupby('abbrev')['conc'].mean().reset_index() if 'abbrev' in loc_data.columns else pd.DataFrame()
        group_data = loc_data.groupby('group')['conc'].mean().reset_index() if 'group' in loc_data.columns else pd.DataFrame()

        pie_html = create_pie_chart_svg(compound_data, group_data)
        if pie_html:
            html_parts.append(pie_html)
            html_parts.append("<hr>")

    # Sample table
    html_parts.append("<table style='width:100%; font-size:12px;'>")
    html_parts.append("<tr><th>Agency</th><th>Year</th><th>Conc</th></tr>")

    for _, row in loc_data.iterrows():
        conc = row['conc']
        conc_text = "B.D." if conc == 0 else f"{conc:.2f}"
        yr_text = int(row['yr']) if pd.notna(row['yr']) else 'N/A'
        agency = row.get('agency', 'N/A')

        html_parts.append(f"<tr><td>{agency}</td><td>{yr_text}</td><td>{conc_text} {unit}</td></tr>")

    html_parts.append("</table>")
    html_parts.append("</div>")

    return ''.join(html_parts)

def aggregate_coords(df, method="Most Recent"):
    """Aggregate data by coordinates"""
    if df.empty:
        return df

    # Clean input data first
    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan)

    # Remove rows with invalid lat/lon or conc
    df = df.dropna(subset=['lat', 'lon'])
    if 'conc' in df.columns:
        df = df[np.isfinite(df['conc'])]

    if method == "Most Recent":
        result = df.sort_values('yr', ascending=False).groupby(['lat', 'lon']).first().reset_index()
    elif method == "Maximum":
        result = df.groupby(['lat', 'lon'])['conc'].max().reset_index()
        result = df.groupby(['lat', 'lon']).first().reset_index().merge(result[['lat', 'lon', 'conc']], on=['lat', 'lon'], suffixes=('_drop', ''))
        result = result[[c for c in result.columns if not c.endswith('_drop')]]
    elif method == "Average":
        result = df.groupby(['lat', 'lon'])['conc'].mean().reset_index()
        result = df.groupby(['lat', 'lon']).first().reset_index().merge(result[['lat', 'lon', 'conc']], on=['lat', 'lon'], suffixes=('_drop', ''))
        result = result[[c for c in result.columns if not c.endswith('_drop')]]
    elif method == "Minimum":
        result = df.groupby(['lat', 'lon'])['conc'].min().reset_index()
        result = df.groupby(['lat', 'lon']).first().reset_index().merge(result[['lat', 'lon', 'conc']], on=['lat', 'lon'], suffixes=('_drop', ''))
        result = result[[c for c in result.columns if not c.endswith('_drop')]]
    else:
        result = df

    # Final cleaning
    result = result.replace([np.inf, -np.inf], np.nan)
    result = result.dropna(subset=['lat', 'lon'])
    if 'conc' in result.columns:
        result = result[np.isfinite(result['conc'])]

    return result

# Load data
all_data = load_data()
spatial_data = load_spatial_data()

agencies = sorted(set().union(*[
    set(all_data[media]['raw']['agency'].unique())
    for media in all_data.keys()
]))

# HEADER
st.markdown("""
<div class="main-header">
    <h1>PFAS in the Delaware River Basin</h1>
    <p>Interactive Data Visualization & Analysis Tool</p>
</div>
""", unsafe_allow_html=True)

# SIDEBAR
with st.sidebar:
    st.markdown("### Data Filters")
    st.markdown("---")

    media_map = {
        "Surface Water": "Water",
        "Groundwater": "GW",
        "Sediment": "Sediment",
        "Tissue": "Tissue"
    }

    dataset = st.selectbox("Media Type", list(media_map.keys()))
    dataset_key = media_map[dataset]

    data_type = st.radio("PFAS Data Type", ["ΣPFAS", "Compounds", "Groups"])
    data_type_key = {"ΣPFAS": "all", "Compounds": "raw", "Groups": "grouped"}[data_type]

    agency_filter = st.selectbox("Agency", ["All"] + agencies)

    current_data = all_data[dataset_key][data_type_key].copy()
    full_data_for_popups = all_data[dataset_key][data_type_key].copy()  # Keep for popups

    if data_type == "Compounds":
        chem_filter = st.selectbox("Compound",
                                   sorted(all_data[dataset_key]['raw']['abbrev'].unique()))
        current_data = current_data[current_data['abbrev'] == chem_filter]
        filter_name = chem_filter
    elif data_type == "Groups":
        group_filter = st.selectbox("Group",
                                    sorted(all_data[dataset_key]['grouped']['group'].unique()))
        current_data = current_data[current_data['group'] == group_filter]
        filter_name = group_filter
    else:
        filter_name = "ΣPFAS"

    if dataset_key == "Tissue" and 'species' in current_data.columns:
        species_list = ["All"] + sorted(current_data['species'].dropna().unique().tolist())
        species_filter = st.selectbox("Species", species_list)
        if species_filter != "All":
            current_data = current_data[current_data['species'] == species_filter]
            full_data_for_popups = full_data_for_popups[full_data_for_popups['species'] == species_filter]

    if 'yr' in current_data.columns and not current_data['yr'].isna().all():
        year_min, year_max = int(current_data['yr'].min()), int(current_data['yr'].max())
        year_range = st.slider("Years", year_min, year_max, (year_min, year_max))
        current_data = current_data[(current_data['yr'] >= year_range[0]) &
                                    (current_data['yr'] <= year_range[1])]
        full_data_for_popups = full_data_for_popups[(full_data_for_popups['yr'] >= year_range[0]) &
                                                      (full_data_for_popups['yr'] <= year_range[1])]

    if agency_filter != "All":
        current_data = current_data[current_data['agency'] == agency_filter]
        full_data_for_popups = full_data_for_popups[full_data_for_popups['agency'] == agency_filter]

    st.markdown("---")
    st.markdown("### Map Options")

    mapped_sample = st.selectbox(
        "Display value",
        ["Most Recent", "Maximum", "Average", "Minimum"],
        help="For locations with multiple years"
    )

    show_huc = st.checkbox("Show HUC12 Averages (based on value displayed)", value=False)
    show_rm = st.checkbox("Show River Miles", value=False)

    st.markdown("---")
    st.markdown("""
    <div class="note-box">
    <strong>📌 Note:</strong><br>
    • Data averaged for multiple samples<br>
    • ΣPFAS = Sum of detected compounds<br>
    • B.D. = Below detection limits
    </div>
    """, unsafe_allow_html=True)

# DATA SUMMARY METRICS
col1, col2, col3, col4 = st.columns(4)

with col1:
    total_samples = len(current_data)
    st.metric("Total Samples", f"{total_samples:,}")

with col2:
    unique_locations = current_data[['lat', 'lon']].drop_duplicates().shape[0]
    st.metric("Locations", f"{unique_locations:,}")

with col3:
    if 'conc' in current_data.columns:
        detected = (current_data['conc'] > 0).sum()
        st.metric("Detections", f"{detected:,}")
    else:
        st.metric("Detections", "N/A")

with col4:
    if 'yr' in current_data.columns and not current_data['yr'].isna().all():
        year_span = f"{int(current_data['yr'].min())}-{int(current_data['yr'].max())}"
        st.metric("Year Range", year_span)
    else:
        st.metric("Year Range", "N/A")

st.markdown("---")

# TABS
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Map", "Estuary Analysis", "EPA Criteria", "History", "About"
])

with tab1:
    st.markdown('<div class="info-card"><h3>Interactive Map</h3><p>Explore PFAS concentrations across the Delaware River Basin</p></div>', unsafe_allow_html=True)

    col_map, col_plots = st.columns([2, 1])

    with col_map:
        # Create map
        m = folium.Map(
            location=[40.61678, -75.22281],
            zoom_start=7,
            tiles='OpenStreetMap'
        )

        # Add DRB boundary
        folium.GeoJson(
            spatial_data['drbbnd'],
            style_function=lambda x: {'color': '#1a3a52', 'weight': 2, 'fillOpacity': 0}
        ).add_to(m)

        # Aggregate data by coordinates
        map_data = aggregate_coords(current_data, mapped_sample)

        # Clean data - remove any rows with invalid values
        if not map_data.empty:
            map_data = map_data.replace([np.inf, -np.inf], np.nan)
            map_data = map_data.dropna(subset=['lat', 'lon', 'conc'])

        unit = "ng/g" if dataset_key in ["Sediment", "Tissue"] else "ng/L"

        # FEATURE 1: HUC12 POLYGONS
        if show_huc and not map_data.empty:
            # Convert to GeoDataFrame
            pts_gdf = gpd.GeoDataFrame(
                map_data,
                geometry=[Point(xy) for xy in zip(map_data['lon'], map_data['lat'])],
                crs="EPSG:4326"
            )

            # Spatial join to assign concentrations to HUC12s
            huc_with_data = spatial_data['huc'].copy()

            # For each HUC, find intersecting points and calculate average
            huc_concentrations = []
            huc_counts = []

            for idx, huc_poly in huc_with_data.iterrows():
                pts_in_huc = pts_gdf[pts_gdf.geometry.within(huc_poly.geometry)]
                if len(pts_in_huc) > 0:
                    avg_conc = pts_in_huc['conc'].mean()
                    huc_concentrations.append(avg_conc)
                    huc_counts.append(len(pts_in_huc))
                else:
                    huc_concentrations.append(np.nan)
                    huc_counts.append(0)

            huc_with_data['conc'] = huc_concentrations
            huc_with_data['n_samples'] = huc_counts

            # Create color scale
            valid_concs = [c for c in huc_concentrations if not np.isnan(c) and np.isfinite(c)]
            if valid_concs:
                colormap = cm.LinearColormap(
                    colors=['#440154', '#31688e', '#35b779', '#fde724'],
                    vmin=min(valid_concs),
                    vmax=max(valid_concs)
                )

                # Add HUC12 polygons
                for idx, row in huc_with_data.iterrows():
                    if not np.isnan(row['conc']) and np.isfinite(row['conc']):
                        folium.GeoJson(
                            row.geometry,
                            style_function=lambda x, conc=row['conc']: {
                                'fillColor': colormap(conc),
                                'color': 'grey',
                                'weight': 1,
                                'fillOpacity': 0.5
                            },
                            popup=folium.Popup(f"""
                                <b>HUC12:</b> {row.get('HU_12_NAME', 'Unknown')}<br>
                                <b>Avg Conc:</b> {row['conc']:.2f} {unit}<br>
                                <b># Locations:</b> {row['n_samples']}<br>
                                <b>Area (sq mi):</b> {row.get('SQMI', 'N/A')}
                            """, max_width=250)
                        ).add_to(m)

        # FEATURE 2: RIVER MILE LABELS
        if show_rm:
            # Add river mile labels for whole numbers
            rm_whole = spatial_data['rm'][spatial_data['rm']['RM'] == spatial_data['rm']['RM'].astype(int)]

            for idx, row in rm_whole.iterrows():
                coords = row.geometry.coords[0]
                folium.Marker(
                    location=[coords[1], coords[0]],
                    icon=folium.DivIcon(html=f"""
                        <div style="
                            font-size: 10px;
                            font-weight: bold;
                            color: #1a3a52;
                            text-shadow: 1px 1px 2px white;
                            background: rgba(255,255,255,0.7);
                            padding: 2px 4px;
                            border-radius: 3px;
                        ">{int(row['RM'])}</div>
                    """)
                ).add_to(m)

        # FEATURE 3 & 4: VARIABLE MARKER SIZES + DETAILED POPUPS
        for idx, row in map_data.iterrows():
            if pd.notna(row['lat']) and pd.notna(row['lon']):
                conc = row['conc']
                color = get_color_for_concentration(conc, dataset_key)
                radius = get_marker_radius(conc, dataset_key)

                # Create detailed popup
                popup_html = create_detailed_popup(
                    row['lat'], row['lon'],
                    full_data_for_popups,
                    data_type_key,
                    dataset_key,
                    filter_name
                )

                folium.CircleMarker(
                    location=[row['lat'], row['lon']],
                    radius=radius,
                    popup=folium.Popup(popup_html, max_width=400),
                    color=color,
                    fill=True,
                    fillColor=color,
                    fillOpacity=0.7,
                    weight=1
                ).add_to(m)

        # Display map
        st_folium(m, width=700, height=600)

    with col_plots:
        # Histogram
        st.markdown("#### Concentration Distribution")

        if not current_data.empty and 'conc' in current_data.columns:
            vals = current_data['conc'].dropna()
            vals_pos = vals[vals > 0]

            if len(vals_pos) > 0:
                fig = go.Figure()
                fig.add_trace(go.Histogram(
                    x=vals_pos,
                    nbinsx=25,
                    marker_color='#2e7aaf',
                    name='Concentration'
                ))
                fig.update_layout(
                    xaxis_title=f"{filter_name} ({unit})",
                    yaxis_title="Count",
                    height=300,
                    margin=dict(l=20, r=20, t=30, b=20),
                    showlegend=False,
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                st.plotly_chart(fig, use_container_width=True)

                # Stats
                st.markdown(f"""
                **Statistics:**
                - N: {len(vals)}
                - B.D.: {(vals == 0).sum()}
                - Mean: {vals_pos.mean():.2f}
                - Median: {vals_pos.median():.2f}
                - Max: {vals_pos.max():.2f}
                """)
            else:
                st.info("No detections in current selection")

        # Year histogram
        st.markdown("#### Sampling Years")

        if not current_data.empty and 'yr' in current_data.columns:
            yrs = current_data['yr'].dropna()

            if len(yrs) > 0:
                fig = go.Figure()
                fig.add_trace(go.Histogram(
                    x=yrs,
                    marker_color='#1a3a52',
                    name='Year'
                ))
                fig.update_layout(
                    xaxis_title="Year",
                    yaxis_title="Count",
                    height=300,
                    margin=dict(l=20, r=20, t=30, b=20),
                    showlegend=False,
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                st.plotly_chart(fig, use_container_width=True)

with tab2:
    st.markdown('<div class="info-card"><h3>Delaware River Estuary Analysis</h3><p>PFAS trends along river miles</p></div>', unsafe_allow_html=True)

    if not current_data.empty:
        # Filter points within Delaware River Bay
        pts_gdf = gpd.GeoDataFrame(
            current_data,
            geometry=[Point(xy) for xy in zip(current_data['lon'], current_data['lat'])],
            crs="EPSG:4326"
        )

        within_bay = pts_gdf[pts_gdf.geometry.within(spatial_data['delrivbay'].unary_union)]

        if not within_bay.empty:
            # Find nearest river mile
            rm_points = spatial_data['rm'].copy()

            nearest_rm = []
            for idx, point in within_bay.iterrows():
                distances = rm_points.distance(point.geometry)
                nearest_idx = distances.idxmin()
                nearest_rm.append(rm_points.loc[nearest_idx, 'RM'])

            within_bay['RM'] = nearest_rm

            # Create plot
            fig = go.Figure()

            # Get unique years and assign colors
            years = sorted(within_bay['yr'].dropna().unique())
            colors = px.colors.sample_colorscale("Spectral", [i/(len(years)-1) if len(years) > 1 else 0 for i in range(len(years))])

            for i, year in enumerate(years):
                year_data = within_bay[within_bay['yr'] == year]
                conc_vals = year_data['conc'].replace(0, np.nan)

                fig.add_trace(go.Scatter(
                    x=year_data['RM'],
                    y=conc_vals,
                    mode='markers',
                    name=str(int(year)),
                    marker=dict(size=10, color=colors[i], opacity=0.7),
                    text=[f"RM: {rm}<br>Conc: {c:.2f}<br>Year: {int(y)}"
                          for rm, c, y in zip(year_data['RM'], year_data['conc'], year_data['yr'])],
                    hovertemplate='%{text}<extra></extra>'
                ))

            fig.update_layout(
                title=f"{filter_name} Along Delaware River Estuary",
                xaxis_title="River Mile",
                yaxis_title=f"{filter_name} ({unit})",
                yaxis_type="log",
                xaxis=dict(autorange="reversed"),
                height=600,
                legend_title="Year",
                hovermode='closest',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )

            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No data points found within the Delaware River Estuary for current selection")
    else:
        st.info("No data available for current selection")

with tab3:
    st.markdown("### EPA Drinking Water Maximum Contaminant Level (MCL) for PFAS")

    mcl_data = {
        "Compound": ["PFOA", "PFOS", "PFHxS", "PFNA", "HFPO-DA (GenX)", "Mixture*"],
        "MCL Goal (ng/L)": ["0", "0", "10", "10", "10", "1 (HI)"],
        "Final MCL (ng/L)": ["4.0", "4.0", "10", "10", "10", "1 (HI)"]
    }
    st.table(pd.DataFrame(mcl_data))
    st.caption("*Mixture of PFHxS, PFNA, HFPO-DA, PFBS (Hazard Index)")

    st.markdown("### EPA Draft Human Health Water Quality Criteria (HHC) for PFAS")

    hhc_data = {
        "Compound": ["PFOA", "PFOS", "PFBS"],
        "Water + Organism HHC (ng/L)": ["0.0009", "0.06", "400"],
        "Organism Only HHC (ng/L)": ["0.0036", "0.07", "500"]
    }
    st.table(pd.DataFrame(hhc_data))

    st.info("Note: These HHC values are draft, non-regulatory criteria. States may adopt their own values.")

with tab4:
    st.markdown("### History of PFAS in the Basin")
    st.markdown("""
    Per- and polyfluoroalkyl substances (PFAS) are a group of more than 13,000 synthetic organic compounds
    characterized by exceptionally strong carbon–fluorine bonds. These bonds make PFAS highly resistant to heat,
    water, and oil, and protect them from natural biological and chemical degradation. As a result, PFAS persist
    for long periods in the environment and are often referred to as 'forever chemicals.'

    The Delaware River Basin holds a unique place in the history of PFAS. It was at DuPont's Chambers Works facility
    in Deepwater, New Jersey, that the first commercially produced PFAS compound—polytetrafluoroethylene (PTFE)—
    was accidentally discovered in 1938, later marketed as Teflon™. This discovery helped establish the region as
    a global center for fluoropolymer research and production.

    This long history of industrial and commercial activity across the Delaware River Basin has left a legacy of
    PFAS pollution that is heavily concentrated in the Lehigh Valley and extends south through the estuary. According
    to the EPA's ECHO database, nearly 3,000 industrial sites in the basin have been identified as potentially involved
    in the production, manufacture, or use of PFAS, with more than 1,500 facilities still active today.
    """)

with tab5:
    st.markdown("### About This Application")
    st.markdown("""
    This interactive application provides basin-wide access to publicly available monitoring data for
    **per- and polyfluoroalkyl substances (PFAS)** in surface water, groundwater, sediment, and biological
    tissue across the Delaware River Basin.

    The app was originally developed by scientists at the **Delaware River Basin Commission (DRBC)** as an
    internal tool to better understand the scope, distribution, and variability of PFAS contamination throughout
    the basin.

    #### Data Sources and Scope
    All data displayed in this application are **publicly available** and were retrieved from:
    - U.S. Environmental Protection Agency's Water Quality Portal
    - U.S. Geological Survey's Water Data API

    The dataset currently includes **over 900 ΣPFAS samples**, most of which are surface and ground water samples
    (>500), in addition to almost 300 tissue samples from 17 different aquatic species.

    #### Contact Information
    - **Jeremy Conkle**, Sr. Chemist/Toxicologist: jeremy.conkle@drbc.gov
    - **Matthew Amato**, Water Resource Scientist: matthew.amato@drbc.gov

    For more information about DRBC's PFAS program, visit:
    [https://www.nj.gov/drbc/programs/quality/pfas.html](https://www.nj.gov/drbc/programs/quality/pfas.html)
    """)

# FOOTER
st.markdown("""
<div class="footer">
    <div class="warning-box">
        ⚠️ <strong>WARNING:</strong> Development of this app is ongoing; information is subject to change.
    </div>
    <p style="margin-top: 1rem; color: #495057;">
        Development funded by EPA's Section 106 Water Pollution Control Grant Program (I-98339317-3)
    </p>
    <p style="font-style: italic; margin-top: 0.5rem; color: #6c757d;">
        App Last Modified: January 16, 2026 | Dataset Last Modified: January 12, 2026
    </p>
</div>
""", unsafe_allow_html=True)
