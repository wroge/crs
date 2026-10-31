package crs

import (
	"math"
)

type LambertConformalConic struct {
	Lonf   float64
	Latf   float64
	Scale  float64
	Eastf  float64
	Northf float64
}

func (cs LambertConformalConic) String() string {
	return build("lambert_conformal_conic").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (lcc LambertConformalConic) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phif := radian(lcc.Latf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	m1 := math.Cos(phif) / math.Sqrt(1-s.E2()*sin2(phif))

	n := math.Sin(phif)
	f := lcc.Scale * m1 / (n * math.Pow(tf, n))
	rf := s.A * f * math.Pow(tf, n)

	phi := radian(lat)

	t := math.Tan(math.Pi/4-phi/2) / math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2)

	r := s.A * f * math.Pow(t, n)
	theta := n * deltaLon(lon, lcc.Lonf)

	east := lcc.Eastf + r*math.Sin(theta)
	north := lcc.Northf + rf - r*math.Cos(theta)

	return east, north, h
}

func (lcc LambertConformalConic) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phif := radian(lcc.Latf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	m1 := math.Cos(phif) / math.Sqrt(1-s.E2()*sin2(phif))

	n := math.Sin(phif)
	f := lcc.Scale * m1 / (n * math.Pow(tf, n))
	rf := s.A * f * math.Pow(tf, n)

	return lccInverse(s, east, north, h, lcc.Eastf, lcc.Northf, n, s.A*f, rf, lcc.Lonf, 0)
}

// lccInverse inverts a Lambert conformal conic.
// af is a*F (Michigan passes the scaled-ellipsoid product).
// theta0 is added to the polar angle (Belgian Lambert).
// n is negative south of the equator; r' takes that sign so the t' power stays real.
func lccInverse(s Spheroid, east, north, h, eastf, northf, n, af, rf, lonf, theta0 float64) (float64, float64, float64) {
	de := east - eastf
	dn := rf - (north - northf)
	ri := math.Hypot(de, dn)
	theta := math.Atan2(de, dn)

	if n < 0 {
		ri = -ri
		theta = math.Atan2(-de, -dn)
	}

	theta += theta0
	ti := math.Pow(ri/af, 1/n)
	phi := math.Pi/2 - 2*math.Atan(ti)
	e := s.E()

	for range 6 {
		esin := e * math.Sin(phi)
		next := math.Pi/2 - 2*math.Atan(ti*math.Pow((1-esin)/(1+esin), e/2))
		if math.Abs(next-phi) < 1e-12 {
			phi = next
			break
		}

		phi = next
	}

	return degree(theta/n + radian(lonf)), degree(phi), h
}
