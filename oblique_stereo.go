package crs

import (
	"math"
)

type ObliqueStereographic struct {
	Lonf, Latf, Scale, Eastf, Northf float64
}

func (cs ObliqueStereographic) String() string {
	return build("oblique_stereographic").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

type obliqueStereoConsts struct {
	R, n, c, chi0, k0, lon0 float64
}

func (cs ObliqueStereographic) consts(s Spheroid) obliqueStereoConsts {
	e := s.E()
	e2 := s.E2()
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	phi0 := radian(cs.Latf)
	lon0 := radian(cs.Lonf)
	sin0 := math.Sin(phi0)
	cos0 := math.Cos(phi0)
	oneMe2Sin2 := 1 - e2*sin0*sin0
	rho0 := s.A * (1 - e2) / math.Pow(oneMe2Sin2, 1.5)
	nu0 := s.A / math.Sqrt(oneMe2Sin2)
	R := math.Sqrt(rho0 * nu0)
	n := math.Sqrt(1 + (e2*math.Pow(cos0, 4))/(1-e2))

	S1 := (1 + sin0) / (1 - sin0)
	S2 := (1 - e*sin0) / (1 + e*sin0)
	w1 := math.Pow(S1*math.Pow(S2, e), n)
	sinChi00 := (w1 - 1) / (w1 + 1)
	c := ((n + sin0) * (1 - sinChi00)) / ((n - sin0) * (1 + sinChi00))
	w2 := c * w1
	chi0 := math.Asin((w2 - 1) / (w2 + 1))

	return obliqueStereoConsts{R: R, n: n, c: c, chi0: chi0, k0: k0, lon0: lon0}
}

func (cs ObliqueStereographic) conformalLat(s Spheroid, phi, c, n float64) float64 {
	e := s.E()
	sinPhi := math.Sin(phi)
	Sa := (1 + sinPhi) / (1 - sinPhi)
	Sb := (1 - e*sinPhi) / (1 + e*sinPhi)
	w := c * math.Pow(Sa*math.Pow(Sb, e), n)

	return math.Asin((w - 1) / (w + 1))
}

func (cs ObliqueStereographic) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	k := cs.consts(s)
	chi := cs.conformalLat(s, radian(lat), k.c, k.n)
	dLam := k.n * (radian(lon) - k.lon0)
	B := 1 + math.Sin(chi)*math.Sin(k.chi0) + math.Cos(chi)*math.Cos(k.chi0)*math.Cos(dLam)
	east := cs.Eastf + 2*k.R*k.k0*math.Cos(chi)*math.Sin(dLam)/B
	north := cs.Northf + 2*k.R*k.k0*(math.Sin(chi)*math.Cos(k.chi0)-math.Cos(chi)*math.Sin(k.chi0)*math.Cos(dLam))/B

	return east, north, h
}

func (cs ObliqueStereographic) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	e := s.E()
	e2 := s.E2()
	k := cs.consts(s)

	de := east - cs.Eastf
	dn := north - cs.Northf

	g := 2 * k.R * k.k0 * math.Tan(math.Pi/4-k.chi0/2)
	hh := 4*k.R*k.k0*math.Tan(k.chi0) + g
	i := math.Atan2(de, hh+dn)
	j := math.Atan2(de, g-dn) - i

	chi := k.chi0 + 2*math.Atan((dn-de*math.Tan(j/2))/(2*k.R*k.k0))
	lon := (j+2*i)/k.n + k.lon0

	psi := math.Log((1+math.Sin(chi))/(k.c*(1-math.Sin(chi)))) / (2 * k.n)
	phi := 2*math.Atan(math.Exp(psi)) - math.Pi/2

	for range 10 {
		sinPhi := math.Sin(phi)
		psiI := math.Log(math.Tan(phi/2+math.Pi/4) * math.Pow((1-e*sinPhi)/(1+e*sinPhi), e/2))
		dPhi := (psiI - psi) * math.Cos(phi) * (1 - e2*sinPhi*sinPhi) / (1 - e2)
		phi -= dPhi
		if math.Abs(dPhi) < 1e-14 {
			break
		}
	}

	return degree(lon), degree(phi), h
}
