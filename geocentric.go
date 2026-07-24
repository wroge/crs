package crs

type Geocentric struct{}

func (g Geocentric) String() string {
	return "geocentric"
}

func (g Geocentric) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return s.GeographicToGeocentric(lon, lat, h)
}

func (g Geocentric) ToGeographic(s Spheroid, x, y, z float64) (float64, float64, float64) {
	return s.GeocentricToGeographic(x, y, z)
}
