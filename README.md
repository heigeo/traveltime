# Travel Time

Description: This tool creates a raster whose cells measure the length of time
in seconds, that it takes water to flow across it and then accumulates
the time from the cell to the outlet of the watershed. water velocity
is calculated as a function of hydraulic radius, Manning's N and slope
Using a DEM as a source, flow direction, flow accumulation and slope (percent)
is calculated. Using flow accumulation, each cell is assigned one of three
flow regimes. Slopes are processed so that there are no cells with zero slopes.

Landcover is used to assign flow N and R values to non channelized cells. Users can
use various landcover sources but these must be matched to the N and R table that
is stored in the /data folder in this distribution. The N and R table must
contain one record for each of the raster valeus in the landcover grid and there
must be the same fields in that table for the program to work. R and N values
are assigned based on landcover and flow regime.

Hydraulic Radius is calucated using (X)a + Y where:
<ul>
<li>X = hydraulic radius factor 1 (User setting, 0.0032 default)  </li>
<li>a = drainage area in square miles (calculated from Flow Accumulation)  </li>
<li>Y = hydraulic radius factor 2 (User setting, 1.7255 default)  </li>
</ul>
	
Velocity is calculated as:  
<ul>
<li>vel(feet/second) = (1.49 * Hydraulic Radius^0.667 * slope^0.5) / Mannings N    </li>
</ul>

Travel Time across individual cells is calculated as:  
<ul>
<li>1 / vel * 0.3048 (if the xy units are meters, otherwise the conversion factor is not used  </li>
</ul> 
