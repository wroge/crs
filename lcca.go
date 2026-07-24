package crs

import (
	"math"
)

type LambertConicNearConformal struct {
	Lonf, Latf, Scale, Eastf, Northf float64
}

func (cs LambertConicNearConformal) String() string {
	return build("lambert_conic_near_conformal").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

type lccaConsts struct {
	l, m0, r0, c, k0 float64
}

func (cs LambertConicNearConformal) consts(s Spheroid) lccaConsts {
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	phi0 := radian(cs.Latf)
	l := math.Sin(phi0)
	m0 := meridianDistance(s, phi0) / s.A
	s2p0 := l * l
	r0 := 1 / (1 - s.E2()*s2p0)
	n0 := math.Sqrt(r0)
	r0 *= (1 - s.E2()) * n0
	tan0 := math.Tan(phi0)

	return lccaConsts{
		l: l, m0: m0, r0: n0 / tan0, c: 1 / (6 * r0 * n0), k0: k0,
	}
}

func lccaFS(S, C float64) float64  { return S * (1 + S*S*C) }
func lccaFSp(S, C float64) float64 { return 1 + 3*S*S*C }

func (cs LambertConicNearConformal) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	q := cs.consts(s)
	S := meridianDistance(s, radian(lat))/s.A - q.m0
	dr := lccaFS(S, q.c)
	r := q.r0 - dr
	lamL := (radian(lon) - radian(cs.Lonf)) * q.l
	east := cs.Eastf + s.A*q.k0*(r*math.Sin(lamL))
	north := cs.Northf + s.A*q.k0*(q.r0-r*math.Cos(lamL))

	return east, north, h
}

func (cs LambertConicNearConformal) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	q := cs.consts(s)
	x := (east - cs.Eastf) / (s.A * q.k0)
	y := (north - cs.Northf) / (s.A * q.k0)
	theta := math.Atan2(x, q.r0-y)
	dr := y - x*math.Tan(0.5*theta)
	S := dr

	for range 10 {
		dif := (lccaFS(S, q.c) - dr) / lccaFSp(S, q.c)
		S -= dif
		if math.Abs(dif) < 1e-12 {
			break
		}
	}

	phi := footpointLatitude(s, s.A*(S+q.m0))
	lam := radian(cs.Lonf) + theta/q.l

	return degree(lam), degree(phi), h
}
