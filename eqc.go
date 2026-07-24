package crs

import (
	"math"
)

type EquidistantCylindrical struct {
	Lonf, Latf, Sp1, Eastf, Northf float64
}

func (cs EquidistantCylindrical) String() string {
	return build("equidistant_cylindrical").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"sp1", cs.Sp1,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs EquidistantCylindrical) rc(s Spheroid) float64 {
	phi1 := radian(cs.Sp1)
	sin1 := math.Sin(phi1)
	nu1 := 1 / math.Sqrt(1-s.E2()*sin1*sin1)

	return nu1 * math.Cos(phi1)
}

func (cs EquidistantCylindrical) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	rc := cs.rc(s)
	m0 := meridianDistance(s, radian(cs.Latf))
	east := cs.Eastf + s.A*rc*(radian(lon)-radian(cs.Lonf))
	north := cs.Northf + meridianDistance(s, radian(lat)) - m0

	return east, north, h
}

func (cs EquidistantCylindrical) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	rc := cs.rc(s)
	m0 := meridianDistance(s, radian(cs.Latf))
	lam := radian(cs.Lonf) + (east-cs.Eastf)/(s.A*rc)
	phi := footpointLatitude(s, north-cs.Northf+m0)

	return degree(lam), degree(phi), h
}
