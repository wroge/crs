package crs

import (
	"math"
)

type AlbersEqualArea struct {
	Lonf   float64
	Latf   float64
	Sp1    float64
	Sp2    float64
	Eastf  float64
	Northf float64
}

func (cs AlbersEqualArea) String() string {
	return build("albers_equal_area").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"sp1", cs.Sp1,
		"sp2", cs.Sp2,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (aea AlbersEqualArea) consts(s Spheroid) (n, c, rf float64) {
	phif := radian(aea.Latf)
	phi1 := radian(aea.Sp1)
	phi2 := radian(aea.Sp2)
	e, e2 := s.E(), s.E2()

	alphaf := authalicQ(math.Sin(phif), e, e2)
	alpha1 := authalicQ(math.Sin(phi1), e, e2)
	alpha2 := authalicQ(math.Sin(phi2), e, e2)

	m1 := math.Cos(phi1) / math.Sqrt(1-e2*sin2(phi1))
	m2 := math.Cos(phi2) / math.Sqrt(1-e2*sin2(phi2))

	if math.Abs(phi1-phi2) < 1e-10 {
		n = math.Sin(phi1)
	} else {
		n = (m1*m1 - m2*m2) / (alpha2 - alpha1)
	}

	c = m1*m1 + n*alpha1
	rf = (s.A * math.Sqrt(c-n*alphaf)) / n

	return n, c, rf
}

func (aea AlbersEqualArea) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	lambdaf := radian(aea.Lonf)
	n, c, rf := aea.consts(s)

	lambda := radian(lon)
	phi := radian(lat)
	alpha := authalicQ(math.Sin(phi), s.E(), s.E2())

	theta := n * (lambda - lambdaf)
	r := (s.A * math.Sqrt(c-n*alpha)) / n

	east := aea.Eastf + r*math.Sin(theta)
	north := aea.Northf + rf - r*math.Cos(theta)

	return east, north, h
}

func (aea AlbersEqualArea) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	lambdaf := radian(aea.Lonf)
	n, c, rf := aea.consts(s)
	qp := authalicQ(1, s.E(), s.E2())

	ri := math.Sqrt(intPow(east-aea.Eastf, 2) + intPow(rf-(north-aea.Northf), 2))
	alphai := (c - (intPow(ri, 2) * intPow(n, 2) / s.A2())) / n
	betai := math.Asin(clamp(alphai/qp, -1, 1))

	var theta float64
	if n > 0 {
		theta = math.Atan2((east - aea.Eastf), (rf - (north - aea.Northf)))
	} else {
		theta = math.Atan2(-(east - aea.Eastf), -(rf - (north - aea.Northf)))
	}

	phi := authalicToGeodetic(betai, s)
	lambda := lambdaf + (theta / n)

	return degree(lambda), degree(phi), h
}
