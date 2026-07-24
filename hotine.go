package crs

import (
	"math"
)

type HotineObliqueMercatorA struct {
	Lonf, Latf, Alpha, Gamma, Scale, Eastf, Northf float64
}

func (cs HotineObliqueMercatorA) String() string {
	return build("hotine_oblique_mercator_a").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"azimuth", cs.Alpha,
		"gamma", cs.Gamma,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

type hotineConsts struct {
	A, B, H, gamma0, lon0, gammaC float64
}

func (cs HotineObliqueMercatorA) consts(s Spheroid) hotineConsts {
	a := s.A
	e := s.E()
	e2 := s.E2()
	kC := cs.Scale
	if kC == 0 {
		kC = 1
	}

	phiC := radian(cs.Latf)
	lambdaC := radian(cs.Lonf)
	alphaC := radian(cs.Alpha)
	gammaC := radian(cs.Gamma)

	sinPhi := math.Sin(phiC)
	cosPhi := math.Cos(phiC)
	sLat := sign(sinPhi)

	B := math.Sqrt(1 + (e2*math.Pow(cosPhi, 4))/(1-e2))
	A := a * B * kC * math.Sqrt(1-e2) / (1 - e2*sinPhi*sinPhi)
	t0 := math.Tan(math.Pi/4-phiC/2) /
		math.Pow((1-e*sinPhi)/(1+e*sinPhi), e/2)
	D := B * math.Sqrt(1-e2) / (cosPhi * math.Sqrt(1-e2*sinPhi*sinPhi))
	D2 := D * D
	if D < 1 {
		D2 = 1
	}
	F := D + math.Sqrt(D2-1)*sLat
	H := F * math.Pow(t0, B)
	G := (F - 1/F) / 2
	gamma0 := math.Asin(clamp(math.Sin(alphaC)/D, -1, 1))
	lon0 := lambdaC - math.Asin(clamp(G*math.Tan(gamma0), -1, 1))/B

	return hotineConsts{A: A, B: B, H: H, gamma0: gamma0, lon0: lon0, gammaC: gammaC}
}

func (cs HotineObliqueMercatorA) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	e := s.E()
	c := cs.consts(s)
	phi := radian(lat)
	lambda := radian(lon)

	t := math.Tan(math.Pi/4-phi/2) /
		math.Pow((1-e*math.Sin(phi))/(1+e*math.Sin(phi)), e/2)
	Q := c.H / math.Pow(t, c.B)
	S := (Q - 1/Q) / 2
	T := (Q + 1/Q) / 2
	V := math.Sin(c.B * (lambda - c.lon0))
	U := (-V*math.Cos(c.gamma0) + S*math.Sin(c.gamma0)) / T

	v := c.A * math.Log((1-U)/(1+U)) / (2 * c.B)
	u := c.A * math.Atan2(
		S*math.Cos(c.gamma0)+V*math.Sin(c.gamma0),
		math.Cos(c.B*(lambda-c.lon0)),
	) / c.B

	east := v*math.Cos(c.gammaC) + u*math.Sin(c.gammaC) + cs.Eastf
	north := u*math.Cos(c.gammaC) - v*math.Sin(c.gammaC) + cs.Northf
	return east, north, h
}

func (cs HotineObliqueMercatorA) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	e2 := s.E2()
	e4 := e2 * e2
	e6 := e4 * e2
	e8 := e4 * e4
	c := cs.consts(s)

	v := (east-cs.Eastf)*math.Cos(c.gammaC) - (north-cs.Northf)*math.Sin(c.gammaC)
	u := (north-cs.Northf)*math.Cos(c.gammaC) + (east-cs.Eastf)*math.Sin(c.gammaC)

	Q := math.Exp(-(c.B * v / c.A))
	S := (Q - 1/Q) / 2
	T := (Q + 1/Q) / 2
	V := math.Sin(c.B * u / c.A)
	U := (V*math.Cos(c.gamma0) + S*math.Sin(c.gamma0)) / T
	t := math.Pow(c.H/math.Sqrt((1+U)/(1-U)), 1/c.B)

	chi := math.Pi/2 - 2*math.Atan(t)
	phi := chi +
		(e2/2+5*e4/24+e6/12+13*e8/360)*math.Sin(2*chi) +
		(7*e4/48+29*e6/240+811*e8/11520)*math.Sin(4*chi) +
		(7*e6/120+81*e8/1120)*math.Sin(6*chi) +
		(4279*e8/161280)*math.Sin(8*chi)

	lon := c.lon0 - math.Atan2(
		S*math.Cos(c.gamma0)-V*math.Sin(c.gamma0),
		math.Cos(c.B*u/c.A),
	)/c.B

	return degree(lon), degree(phi), h
}
