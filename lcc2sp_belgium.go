package crs

import (
	"math"
)

const belgiumLCCAngle = 29.2985 * math.Pi / (180 * 3600)

type LambertConformalConic2SPBelgium struct {
	Lonf, Latf, Sp1, Sp2, Eastf, Northf float64
}

func (cs LambertConformalConic2SPBelgium) String() string {
	return build("lambert_conformal_conic_2sp_belgium").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"sp1", cs.Sp1,
		"sp2", cs.Sp2,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LambertConformalConic2SPBelgium) consts(s Spheroid) (n, f, rf float64) {
	phif := radian(cs.Latf)
	phi1 := radian(cs.Sp1)
	phi2 := radian(cs.Sp2)
	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	t1 := math.Tan(math.Pi/4-phi1/2) / math.Pow((1-s.E()*math.Sin(phi1))/(1+s.E()*math.Sin(phi1)), s.E()/2)
	t2 := math.Tan(math.Pi/4-phi2/2) / math.Pow((1-s.E()*math.Sin(phi2))/(1+s.E()*math.Sin(phi2)), s.E()/2)
	m1 := math.Cos(phi1) / math.Sqrt(1-s.E2()*sin2(phi1))
	m2 := math.Cos(phi2) / math.Sqrt(1-s.E2()*sin2(phi2))

	if math.Abs(phi1-phi2) < 1e-14 {
		n = math.Sin(phi1)
	} else {
		n = math.Log(m1/m2) / math.Log(t1/t2)
	}

	f = m1 / (n * math.Pow(t1, n))
	rf = s.A * f * math.Pow(tf, n)

	return n, f, rf
}

func (cs LambertConformalConic2SPBelgium) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	n, f, rf := cs.consts(s)
	phi := radian(lat)
	t := math.Tan(math.Pi/4-phi/2) / math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2)
	r := s.A * f * math.Pow(t, n)
	theta := n*(radian(lon)-radian(cs.Lonf)) - belgiumLCCAngle
	east := cs.Eastf + r*math.Sin(theta)
	north := cs.Northf + rf - r*math.Cos(theta)

	return east, north, h
}

func (cs LambertConformalConic2SPBelgium) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	n, f, rf := cs.consts(s)
	ri := math.Hypot(east-cs.Eastf, rf-(north-cs.Northf))
	ti := math.Pow(ri/(s.A*f), 1/n)
	theta := math.Atan2(east-cs.Eastf, rf-(north-cs.Northf)) + belgiumLCCAngle

	phi := math.Pi/2 - 2*math.Atan(ti)

	for range 6 {
		next := math.Pi/2 - 2*math.Atan(ti*math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2))
		if math.Abs(next-phi) < 1e-12 {
			phi = next
			break
		}
		phi = next
	}

	return degree(theta/n + radian(cs.Lonf)), degree(phi), h
}
