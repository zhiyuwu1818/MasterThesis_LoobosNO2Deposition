import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from pyproj import Transformer

def export_footprint_shapefile(FFP, year_month, rs, lat, lon, output_folder = '../Data/Footprint'):
    """
    Export the footprint contour as a shapefile.

    Parameters:
    FFP : dict
        Footprint data containing 'x', 'y', and 'FFP' arrays.
    period_label : str
        Descriptive label for the period being exported. Can be any string,
        e.g., '2023' (annual), '2023-07' (monthly), 'spring_2023' (seasonal).
    rs : float
        The footprint contour level to export (e.g., 0.8 for 80%).
    lat : float
        Latitude of the flux tower.
    lon : float
        Longitude of the flux tower.
    output_folder : str
        Folder to save the output shapefile.
    """
    # Convert lat and lon to RD new coordinates
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:28992", always_xy=True)
    tower_x, tower_y = transformer.transform(lon, lat)    
    # Convert FFP x and y to absolute coordinates
    X = np.array(FFP['xr'][0]) + tower_x
    Y = np.array(FFP['yr'][0]) + tower_y
    # Create the footprint contour polygon
    footprint_polygon = Polygon(zip(X, Y))
    # Create a GeoDataFrame, and then transform to EPSG:32631 for saving
    gdf = gpd.GeoDataFrame({'contour level': [rs]}, crs="EPSG:28992", geometry=[footprint_polygon])
    gdf = gdf.to_crs("EPSG:32631")

    # Save to shapefile
    output_path = f"{output_folder}/footprint_{year_month}_{int(rs*100)}pct_S2CRS.shp"
    gdf.to_file(output_path)
    print(f"Footprint contour shapefile saved to: {output_path}")