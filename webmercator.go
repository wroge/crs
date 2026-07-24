package crs

import (
	"math"
)

type WebMercator struct {
	Lonf   float64
	Latf   float64
	Scale  float64
	Eastf  float64
	Northf float64
}

func (cs WebMercator) String() string {
	return build("web_mercator").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs WebMercator) k() float64 {
	if cs.Scale == 0 {
		return 1
	}

	return cs.Scale
}

func mercatorM(phi float64) float64 {
	return math.Log(math.Tan(math.Pi/4 + phi/2))
}

func (cs WebMercator) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	lambda := radian(lon)
	phi := radian(lat)
	lambda0 := radian(cs.Lonf)
	phi0 := radian(cs.Latf)
	k := cs.k()
	rk := k * s.A

	east := cs.Eastf + rk*(lambda-lambda0)
	north := cs.Northf + rk*(mercatorM(phi)-mercatorM(phi0))

	return east, north, h
}

func (cs WebMercator) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	lambda0 := radian(cs.Lonf)
	phi0 := radian(cs.Latf)
	k := cs.k()
	rk := k * s.A

	psi := (north-cs.Northf)/rk + mercatorM(phi0)
	phi := 2*math.Atan(math.Exp(psi)) - math.Pi/2
	lambda := lambda0 + (east-cs.Eastf)/rk

	return degree(lambda), degree(phi), h
}
