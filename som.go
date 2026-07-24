package crs

import (
	"math"
)

// SwissObliqueMercator is EPSG method 9814 (PROJ +proj=somerc).
// Alpha/Gamma are retained for EPSG parameter round-trip; the classic Swiss
// formulas assume azimuth = rectified grid angle = 90°.
type SwissObliqueMercator struct {
	Lonf   float64
	Latf   float64
	Scale  float64
	Eastf  float64
	Northf float64
	Alpha  float64
	Gamma  float64
}

func (cs SwissObliqueMercator) String() string {
	return build("swiss_oblique_mercator").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
		"alpha", cs.Alpha,
		"gamma", cs.Gamma,
	).String()
}

type somercConsts struct {
	c, K, kR, cosp0, sinp0, hlfE, lon0 float64
}

func (som SwissObliqueMercator) consts(s Spheroid) somercConsts {
	e := s.E()
	e2 := s.E2()
	k0 := som.Scale
	if k0 == 0 {
		k0 = 1
	}

	phi0 := radian(som.Latf)
	cos0 := math.Cos(phi0)
	cos02 := cos0 * cos0
	c := math.Sqrt(1 + e2*cos02*cos02/(1-e2))
	sp := math.Sin(phi0)
	sinp0 := sp / c
	phip0 := math.Asin(clamp(sinp0, -1, 1))
	cosp0 := math.Cos(phip0)
	esp := e * sp
	K := math.Log(math.Tan(math.Pi/4+phip0/2)) -
		c*(math.Log(math.Tan(math.Pi/4+phi0/2))-0.5*e*math.Log((1+esp)/(1-esp)))
	kR := s.A * k0 * math.Sqrt(1-e2) / (1 - esp*esp)

	return somercConsts{
		c: c, K: K, kR: kR, cosp0: cosp0, sinp0: sinp0, hlfE: 0.5 * e, lon0: radian(som.Lonf),
	}
}

func (som SwissObliqueMercator) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	q := som.consts(s)
	phi := radian(lat)
	lam := radian(lon) - q.lon0

	sp := s.E() * math.Sin(phi)
	phip := 2*math.Atan(math.Exp(
		q.c*(math.Log(math.Tan(math.Pi/4+phi/2))-q.hlfE*math.Log((1+sp)/(1-sp)))+q.K,
	)) - math.Pi/2
	lamp := q.c * lam
	cp := math.Cos(phip)
	phipp := math.Asin(clamp(q.cosp0*math.Sin(phip)-q.sinp0*cp*math.Cos(lamp), -1, 1))
	lampp := math.Asin(clamp(cp*math.Sin(lamp)/math.Cos(phipp), -1, 1))

	east := q.kR*lampp + som.Eastf
	north := q.kR*math.Log(math.Tan(math.Pi/4+phipp/2)) + som.Northf

	return east, north, h
}

func (som SwissObliqueMercator) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	q := som.consts(s)
	e := s.E()
	roneEs := 1 / (1 - s.E2())

	phipp := 2 * (math.Atan(math.Exp((north-som.Northf)/q.kR)) - math.Pi/4)
	lampp := (east - som.Eastf) / q.kR
	cp := math.Cos(phipp)
	phip := math.Asin(clamp(q.cosp0*math.Sin(phipp)+q.sinp0*cp*math.Cos(lampp), -1, 1))
	lamp := math.Asin(clamp(cp*math.Sin(lampp)/math.Cos(phip), -1, 1))

	con := (q.K - math.Log(math.Tan(math.Pi/4+phip/2))) / q.c
	phi := phip
	for range 6 {
		esp := e * math.Sin(phi)
		delp := (con + math.Log(math.Tan(math.Pi/4+phi/2)) -
			q.hlfE*math.Log((1+esp)/(1-esp))) *
			(1 - esp*esp) * math.Cos(phi) * roneEs
		phi -= delp
		if math.Abs(delp) < 1e-10 {
			break
		}
	}

	return degree(q.lon0 + lamp/q.c), degree(phi), h
}
