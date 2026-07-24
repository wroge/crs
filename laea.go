package crs

import (
	"math"
)

type LambertAzimuthalEqualArea struct {
	Lonf   float64
	Latf   float64
	Eastf  float64
	Northf float64
}

func (cs LambertAzimuthalEqualArea) String() string {
	return build("lambert_azimuthal_equal_area").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (laea LambertAzimuthalEqualArea) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phi0 := radian(laea.Latf)
	lambda0 := radian(laea.Lonf)
	e, e2 := s.E(), s.E2()

	q0 := authalicQ(math.Sin(phi0), e, e2)
	qp := authalicQ(1, e, e2)

	beta0 := math.Asin(q0 / qp)
	rq := s.A * math.Sqrt(qp/2)
	g := s.A * (math.Cos(phi0) / math.Sqrt(1-e2*sin2(phi0))) / (rq * math.Cos(beta0))

	phi := radian(lat)
	lambda := radian(lon)

	q := authalicQ(math.Sin(phi), e, e2)
	beta := math.Asin(q / qp)
	b := rq * math.Sqrt(2/(1+math.Sin(beta0)*math.Sin(beta)+(math.Cos(beta0)*math.Cos(beta)*math.Cos(lambda-lambda0))))

	east := laea.Eastf + ((b * g) * (math.Cos(beta) * math.Sin(lambda-lambda0)))
	north := laea.Northf + (b/g)*((math.Cos(beta0)*math.Sin(beta))-(math.Sin(beta0)*math.Cos(beta)*math.Cos(lambda-lambda0)))

	return east, north, h
}

func (laea LambertAzimuthalEqualArea) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phi0 := radian(laea.Latf)
	lambda0 := radian(laea.Lonf)
	e, e2 := s.E(), s.E2()

	q0 := authalicQ(math.Sin(phi0), e, e2)
	qp := authalicQ(1, e, e2)

	beta0 := math.Asin(q0 / qp)
	rq := s.A * math.Sqrt(qp/2)
	g := s.A * (math.Cos(phi0) / math.Sqrt(1-e2*sin2(phi0))) / (rq * math.Cos(beta0))

	rho := math.Sqrt(intPow((east-laea.Eastf)/g, 2) + intPow(g*(north-laea.Northf), 2))

	if rho < 1e-14 {
		return laea.Lonf, laea.Latf, h
	}

	c := 2 * math.Asin(rho/(2*rq))
	betai := math.Asin((math.Cos(c) * math.Sin(beta0)) + ((g * (north - laea.Northf) * math.Sin(c) * math.Cos(beta0)) / rho))

	phi := authalicToGeodetic(betai, s)
	lambda := lambda0 + math.Atan2((east-laea.Eastf)*math.Sin(c), (g*rho*math.Cos(beta0)*math.Cos(c)-intPow(g, 2)*(north-laea.Northf)*math.Sin(beta0)*math.Sin(c)))

	return degree(lambda), degree(phi), h
}
