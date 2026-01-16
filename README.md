# Delaware River Basin PFAS Data Visualization App
## Python/Streamlit Version 🐍

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-app.streamlit.app)

![PFAS App Screenshot](https://via.placeholder.com/800x400?text=PFAS+Data+Visualization+App)

## Overview

This is a complete Python/Streamlit port of the Delaware River Basin PFAS (Per- and polyfluoroalkyl substances) Data Visualization Application. The app provides interactive access to publicly available monitoring data for PFAS across surface water, groundwater, sediment, and biological tissue in the Delaware River Basin.

## ✨ Features

### 🗺️ Interactive Map
- Folium-based interactive mapping
- Color-coded concentration markers
- Clickable popups with sample details
- Delaware River Basin boundary overlay
- Multiple data aggregation methods (Most Recent, Maximum, Average, Minimum)

### 📊 Data Visualizations
- Real-time concentration histograms with Plotly
- Year distribution analysis
- River mile trend analysis (estuary-specific)
- Summary statistics dashboard

### 🎛️ Advanced Filtering
- **Media Types**: Surface Water, Groundwater, Sediment, Tissue
- **Data Types**: ΣPFAS (sum), Individual Compounds, Chemical Groups
- **Agency Filter**: DRBC, USGS, USEPA, NJDEP, Other
- **Year Range Slider**: Filter by sampling year
- **Species Filter**: For tissue data (17 aquatic species)
- **Compound/Group Selection**: Dynamic based on data type

### 📋 Additional Resources
- EPA Drinking Water MCL criteria tables
- EPA Draft Human Health Water Quality Criteria
- PFAS history in the Delaware River Basin
- Complete documentation and data sources

### 🎨 Modern UI
- Beautiful gradient design theme
- Responsive layout
- Professional styling
- Real-time metric cards
- Interactive Plotly charts (zoom, pan, export)

## 🚀 Quick Start

### Deploy to Streamlit Cloud (5 minutes - FREE!)

1. **Fork this repository** or use it directly

2. **Go to [Streamlit Cloud](https://share.streamlit.io)**

3. **Create new app**:
   - Repository: `acktek/DRB-PFAS-Python`
   - Branch: `main`
   - Main file path: `app.py`

4. **Click Deploy** - Your app will be live in 2-3 minutes!

### Run Locally

```bash
# Clone the repository
git clone https://github.com/acktek/DRB-PFAS-Python.git
cd DRB-PFAS-Python

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

## 📦 Dependencies

- Python 3.11+
- streamlit
- pandas
- geopandas
- folium
- plotly
- shapely
- numpy

See `requirements.txt` for complete list.

## 📁 Data Sources

All data displayed in this application are **publicly available** from:

- **U.S. EPA Water Quality Portal**
- **U.S. Geological Survey Water Data API**

The dataset includes:
- **900+ ΣPFAS samples**
- **500+ water samples** (surface water and groundwater)
- **300+ tissue samples** from 17 different aquatic species
- Spatial data for the Delaware River Basin

## 🏗️ Project Structure

```
DRB-PFAS-Python/
├── app.py                              # Main Streamlit application
├── requirements.txt                     # Python dependencies
├── Dockerfile                          # Docker configuration
├── DEPLOY.md                           # Deployment guide
├── PFAS_water_data_DRB_20260112.csv   # Surface water data
├── PFAS_ground_water_data_DRB_20260112.csv  # Groundwater data
├── PFAS_sediment_data_DRB_20260112.csv      # Sediment data
├── PFAS_tissue_data_DRB_20260112.csv        # Tissue data
├── *.shp, *.shx, *.dbf, *.prj        # Shapefiles (spatial data)
└── www/                                # Static assets (logos, images)
```

## 🌐 Deployment Options

### Streamlit Cloud (Recommended - FREE)
- ✅ Completely free, unlimited hosting
- ✅ Auto-deploys from GitHub
- ✅ No server management
- ✅ Built-in SSL/HTTPS
- ✅ Custom subdomain

### Render
- ✅ Free tier available
- ✅ Docker support
- ⚠️ May sleep after inactivity

### Hugging Face Spaces
- ✅ Free, unlimited
- ✅ Great for data/ML apps
- ✅ Easy deployment

## 🔄 Comparison with R Shiny Version

| Feature | R Shiny | Python/Streamlit |
|---------|---------|------------------|
| Interactive Map | ✅ Leaflet | ✅ Folium |
| Histograms | ✅ Base R plots | ✅ Interactive Plotly |
| River Mile Analysis | ✅ Base R plots | ✅ Interactive Plotly |
| Filters & Controls | ✅ | ✅ |
| EPA Criteria | ✅ | ✅ |
| Beautiful UI | ✅ | ✅ Enhanced gradients |
| Deployment | Render (sleeps) | Streamlit (always on!) |
| **Cost** | **FREE** | **FREE** |

**Bonus in Python version:**
- Interactive zoom/pan on charts
- Faster initial load
- Easier deployment
- Modern gradient UI
- Export chart capabilities

## 📝 Notes

- Data are averaged when multiple samples exist at the same location within the same calendar year
- ΣPFAS = Sum of detected PFAS compounds in a sample
- B.D. = Below analytical detection limits (varies by sample)
- Not all samples were analyzed for the same compounds

## 👥 Contact

This application was developed by scientists at the **Delaware River Basin Commission (DRBC)**.

- **Jeremy Conkle**, Sr. Chemist/Toxicologist: [jeremy.conkle@drbc.gov](mailto:jeremy.conkle@drbc.gov)
- **Matthew Amato**, Water Resource Scientist: [matthew.amato@drbc.gov](mailto:matthew.amato@drbc.gov)

For more information about DRBC's PFAS program:
[https://www.nj.gov/drbc/programs/quality/pfas.html](https://www.nj.gov/drbc/programs/quality/pfas.html)

## 📄 License

Development funded by EPA's Section 106 Water Pollution Control Grant Program (I-98339317-3)

## ⚠️ Disclaimer

**WARNING**: Development of this app is ongoing; information is subject to change.

**Last Updated**: January 16, 2026
**Dataset Last Modified**: January 12, 2026

---

## 🤝 Contributing

Issues and pull requests are welcome! Please feel free to contribute to improving this application.

## 🌟 Acknowledgments

- Delaware River Basin Commission (DRBC)
- U.S. Environmental Protection Agency (EPA)
- U.S. Geological Survey (USGS)

---

Made with ❤️ using Python and Streamlit
