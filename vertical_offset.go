package crs

import (
	"math"
)

type VerticalOffset struct {
	Dh float64
}

func (o VerticalOffset) String() string {
	return build("").addAll(
		"operation", "vertical_offset",
		"dh", o.Dh,
	).String()
}

func (o VerticalOffset) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return lon, lat, h + o.Dh, nil
}

func (o VerticalOffset) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return lon0, lat0, h0 - o.Dh, nil
}

type VerticalOffsetAndSlope struct {
	Lat0, Lon0         float64
	Dh                 float64
	SlopeLat, SlopeLon float64 // arc-seconds
}

func (o VerticalOffsetAndSlope) String() string {
	return build("").addAll(
		"operation", "vertical_offset_and_slope",
		"lat0", o.Lat0,
		"lon0", o.Lon0,
		"dh", o.Dh,
		"slope_lat", o.SlopeLat,
		"slope_lon", o.SlopeLon,
	).String()
}

func (o VerticalOffsetAndSlope) delta(s Spheroid, lon, lat float64) float64 {
	phi0 := radian(o.Lat0)
	lam0 := radian(o.Lon0)
	phi := radian(lat)
	lam := radian(lon)

	sinPhi0 := math.Sin(phi0)
	e2 := s.E2()
	oneMe2Sin2 := 1 - e2*sinPhi0*sinPhi0
	rho0 := s.A * (1 - e2) / (oneMe2Sin2 * math.Sqrt(oneMe2Sin2))
	nu0 := s.A / math.Sqrt(oneMe2Sin2)

	slopeLat := radian(o.SlopeLat / 3600)
	slopeLon := radian(o.SlopeLon / 3600)

	return o.Dh +
		slopeLat*rho0*(phi-phi0) +
		slopeLon*nu0*(lam-lam0)*math.Cos(phi)
}

// ToTarget implements [Operation].
func (o VerticalOffsetAndSlope) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return lon, lat, h + o.delta(source, lon, lat), nil
}

// FromTarget implements [Operation].
func (o VerticalOffsetAndSlope) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return lon0, lat0, h0 - o.delta(source, lon0, lat0), nil
}
