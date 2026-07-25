package crs

import (
	"math"
)

type HyperbolicCassiniSoldner struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs HyperbolicCassiniSoldner) String() string {
	return build("hyperbolic_cassini_soldner").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs HyperbolicCassiniSoldner) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	base := CassiniSoldner(cs)

	east, north, h := base.FromGeographic(s, lon, lat, h)
	phi := radian(lat)
	sinPhi := math.Sin(phi)
	e2 := s.E2()
	oneMe2Sin2 := 1 - e2*sinPhi*sinPhi
	nu := s.A / math.Sqrt(oneMe2Sin2)
	rho := s.A * (1 - e2) / math.Pow(oneMe2Sin2, 1.5)
	x := north - cs.Northf
	north = cs.Northf + x - x*x*x/(6*rho*nu)

	return east, north, h
}

func (cs HyperbolicCassiniSoldner) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phi0 := radian(cs.Latf)
	target := north - cs.Northf
	x := target
	e2 := s.E2()

	for range 10 {
		m1 := meridianDistance(s, phi0) + x
		phi1 := footpointLatitude(s, m1)
		sinPhi1 := math.Sin(phi1)
		oneMe2Sin2 := 1 - e2*sinPhi1*sinPhi1
		nu := s.A / math.Sqrt(oneMe2Sin2)
		rho := s.A * (1 - e2) / math.Pow(oneMe2Sin2, 1.5)
		f := x - x*x*x/(6*rho*nu) - target
		df := 1 - x*x/(2*rho*nu)
		dx := f / df
		x -= dx
		if math.Abs(dx) < 1e-12 {
			break
		}
	}

	base := CassiniSoldner(cs)

	return base.ToGeographic(s, east, cs.Northf+x, h)
}
