package crs

import (
	"math"
)

// MercatorA is EPSG method 9804 (PROJ +proj=merc with k / 1SP).
// Latf is unused in the formulas (natural origin is on the equator).
type MercatorA struct {
	Lonf, Latf, Scale, Eastf, Northf float64
}

func (cs MercatorA) String() string {
	return build("mercator_a").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs MercatorA) k0() float64 {
	if cs.Scale == 0 {
		return 1
	}
	return cs.Scale
}

func (cs MercatorA) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return mercatorFromGeographic(s, lon, lat, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0())
}

func (cs MercatorA) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	return mercatorToGeographic(s, east, north, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0())
}

// MercatorB is EPSG method 9805 (PROJ +proj=merc with lat_ts / 2SP).
// Sp1 is the latitude of the first standard parallel; k0 is derived from it.
type MercatorB struct {
	Lonf, Sp1, Eastf, Northf float64
}

func (cs MercatorB) String() string {
	return build("mercator_b").addAll(
		"lonf", cs.Lonf,
		"sp1", cs.Sp1,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs MercatorB) k0(s Spheroid) float64 {
	phi1 := radian(cs.Sp1)
	sin1 := math.Sin(phi1)
	return math.Cos(phi1) / math.Sqrt(1-s.E2()*sin1*sin1)
}

func (cs MercatorB) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return mercatorFromGeographic(s, lon, lat, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(s))
}

func (cs MercatorB) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	return mercatorToGeographic(s, east, north, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(s))
}

func mercatorFromGeographic(s Spheroid, lon, lat, h, lonf, eastf, northf, k0 float64) (float64, float64, float64) {
	e := s.E()
	phi := radian(lat)
	sinPhi := math.Sin(phi)

	east := eastf + s.A*k0*(radian(lon)-radian(lonf))
	north := northf + s.A*k0*math.Log(
		math.Tan(math.Pi/4+phi/2)*math.Pow((1-e*sinPhi)/(1+e*sinPhi), e/2),
	)

	return east, north, h
}

func mercatorToGeographic(s Spheroid, east, north, h, lonf, eastf, northf, k0 float64) (float64, float64, float64) {
	e2 := s.E2()
	e4 := e2 * e2
	e6 := e4 * e2
	e8 := e4 * e4

	t := math.Exp((northf - north) / (s.A * k0))
	chi := math.Pi/2 - 2*math.Atan(t)

	phi := chi +
		(e2/2+5*e4/24+e6/12+13*e8/360)*math.Sin(2*chi) +
		(7*e4/48+29*e6/240+811*e8/11520)*math.Sin(4*chi) +
		(7*e6/120+81*e8/1120)*math.Sin(6*chi) +
		(4279*e8/161280)*math.Sin(8*chi)
	lambda := radian(lonf) + (east-eastf)/(s.A*k0)

	return degree(lambda), degree(phi), h
}
