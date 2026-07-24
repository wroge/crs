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
	lambdaf := radian(lcc.Lonf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	m1 := math.Cos(phif) / math.Sqrt(1-s.E2()*sin2(phif))

	n := math.Sin(phif)
	f := lcc.Scale * m1 / (n * math.Pow(tf, n))
	rf := s.A * f * math.Pow(tf, n)

	phi := radian(lat)
	lambda := radian(lon)

	t := math.Tan(math.Pi/4-phi/2) / math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2)

	r := s.A * f * math.Pow(t, n)
	theta := n * (lambda - lambdaf)

	east := lcc.Eastf + r*math.Sin(theta)
	north := lcc.Northf + rf - r*math.Cos(theta)

	return east, north, h
}

func (lcc LambertConformalConic) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phif := radian(lcc.Latf)
	lambdaf := radian(lcc.Lonf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	m1 := math.Cos(phif) / math.Sqrt(1-s.E2()*sin2(phif))

	n := math.Sin(phif)
	f := lcc.Scale * m1 / (n * math.Pow(tf, n))
	rf := s.A * f * math.Pow(tf, n)
	ri := math.Hypot(east-lcc.Eastf, rf-(north-lcc.Northf))

	ti := math.Pow(ri/(s.A*f), 1/n)

	theta := math.Atan2(east-lcc.Eastf, rf-(north-lcc.Northf))

	phi := math.Pi/2 - 2*math.Atan(ti)

	for range 6 {
		next := math.Pi/2 - 2*math.Atan(ti*math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2))
		if math.Abs(next-phi) < 1e-12 {
			phi = next
			break
		}

		phi = next
	}

	lambda := theta/n + lambdaf

	return degree(lambda), degree(phi), h
}
