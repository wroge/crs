package crs

type Geographic struct{}

func (g Geographic) String() string {
	return "geographic"
}

func (g Geographic) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return lon, lat, h
}

func (g Geographic) ToGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return lon, lat, h
}
