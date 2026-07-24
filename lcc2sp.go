package crs

import (
	"math"
)

type LambertConformalConic2SP struct {
	Lonf   float64
	Latf   float64
	Sp1    float64
	Sp2    float64
	Eastf  float64
	Northf float64
}

func (cs LambertConformalConic2SP) String() string {
	return build("lambert_conformal_conic_2sp").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"sp1", cs.Sp1,
		"sp2", cs.Sp2,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (lcc LambertConformalConic2SP) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phif := radian(lcc.Latf)
	phi1 := radian(lcc.Sp1)
	phi2 := radian(lcc.Sp2)
	lambdaf := radian(lcc.Lonf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	t1 := math.Tan(math.Pi/4-phi1/2) / math.Pow((1-s.E()*math.Sin(phi1))/(1+s.E()*math.Sin(phi1)), s.E()/2)
	t2 := math.Tan(math.Pi/4-phi2/2) / math.Pow((1-s.E()*math.Sin(phi2))/(1+s.E()*math.Sin(phi2)), s.E()/2)

	m1 := math.Cos(phi1) / math.Sqrt(1-s.E2()*sin2(phi1))
	m2 := math.Cos(phi2) / math.Sqrt(1-s.E2()*sin2(phi2))

	var n float64
	if math.Abs(phi1-phi2) < 1e-14 {
		n = math.Sin(phi1)
	} else {
		n = math.Log(m1/m2) / math.Log(t1/t2)
	}

	f := m1 / (n * math.Pow(t1, n))
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

func (lcc LambertConformalConic2SP) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phif := radian(lcc.Latf)
	phi1 := radian(lcc.Sp1)
	phi2 := radian(lcc.Sp2)
	lambdaf := radian(lcc.Lonf)

	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	t1 := math.Tan(math.Pi/4-phi1/2) / math.Pow((1-s.E()*math.Sin(phi1))/(1+s.E()*math.Sin(phi1)), s.E()/2)
	t2 := math.Tan(math.Pi/4-phi2/2) / math.Pow((1-s.E()*math.Sin(phi2))/(1+s.E()*math.Sin(phi2)), s.E()/2)

	m1 := math.Cos(phi1) / math.Sqrt(1-s.E2()*sin2(phi1))
	m2 := math.Cos(phi2) / math.Sqrt(1-s.E2()*sin2(phi2))

	var n float64
	if math.Abs(phi1-phi2) < 1e-14 {
		n = math.Sin(phi1)
	} else {
		n = math.Log(m1/m2) / math.Log(t1/t2)
	}

	f := m1 / (n * math.Pow(t1, n))
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
