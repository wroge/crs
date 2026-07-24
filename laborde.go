package crs

import (
	"math"
)

type LabordeObliqueMercator struct {
	Lonf, Latf, Alpha, Scale, Eastf, Northf float64
}

func (cs LabordeObliqueMercator) String() string {
	return build("laborde_oblique_mercator").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"azimuth", cs.Alpha,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

type labordeConsts struct {
	kRg, p0s, A, C, Ca, Cb, Cc, Cd float64
}

func (cs LabordeObliqueMercator) consts(s Spheroid) labordeConsts {
	k0 := cs.Scale
	if k0 == 0 {
		k0 = 1
	}
	phi0 := radian(cs.Latf)
	az := radian(cs.Alpha)
	e := s.E()
	e2 := s.E2()
	sinp := math.Sin(phi0)
	t := 1 - e2*sinp*sinp
	N := 1 / math.Sqrt(t)
	R := (1 - e2) * N / t

	kRg := k0 * math.Sqrt(N*R)
	p0s := math.Atan(math.Sqrt(R/N) * math.Tan(phi0))
	A := sinp / math.Sin(p0s)
	te := e * sinp
	C := 0.5*e*A*math.Log((1+te)/(1-te)) - A*math.Log(math.Tan(math.Pi/4+0.5*phi0)) +
		math.Log(math.Tan(math.Pi/4+0.5*p0s))
	tAz := az + az
	Cb := 1 / (12 * kRg * kRg)
	Ca := (1 - math.Cos(tAz)) * Cb
	Cb *= math.Sin(tAz)
	Cc := 3 * (Ca*Ca - Cb*Cb)
	Cd := 6 * Ca * Cb

	return labordeConsts{kRg: kRg, p0s: p0s, A: A, C: C, Ca: Ca, Cb: Cb, Cc: Cc, Cd: Cd}
}

func (cs LabordeObliqueMercator) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	q := cs.consts(s)
	e := s.E()
	phi := radian(lat)
	lam := radian(lon - cs.Lonf)

	V1 := q.A * math.Log(math.Tan(math.Pi/4+0.5*phi))
	t := e * math.Sin(phi)
	V2 := 0.5 * e * q.A * math.Log((1+t)/(1-t))
	ps := 2 * (math.Atan(math.Exp(V1-V2+q.C)) - math.Pi/4)
	I1 := ps - q.p0s
	cosps := math.Cos(ps)
	sinps := math.Sin(ps)
	cosps2 := cosps * cosps
	sinps2 := sinps * sinps
	I4 := q.A * cosps
	I2 := 0.5 * q.A * I4 * sinps
	I3 := I2 * q.A * q.A * (5*cosps2 - sinps2) / 12
	I6 := I4 * q.A * q.A
	I5 := I6 * (cosps2 - sinps2) / 6
	I6 *= q.A * q.A * (5*cosps2*cosps2 + sinps2*(sinps2-18*cosps2)) / 120
	t2 := lam * lam
	x := q.kRg * lam * (I4 + t2*(I5+t2*I6))
	y := q.kRg * (I1 + t2*(I2+t2*I3))
	x2 := x * x
	y2 := y * y
	V1 = 3*x*y2 - x*x2
	V2 = y*y2 - 3*x2*y
	x += q.Ca*V1 + q.Cb*V2
	y += q.Ca*V2 - q.Cb*V1

	return cs.Eastf + s.A*x, cs.Northf + s.A*y, h
}

func (cs LabordeObliqueMercator) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	q := cs.consts(s)
	e := s.E()
	e2 := s.E2()
	phi0 := radian(cs.Latf)
	x := (east - cs.Eastf) / s.A
	y := (north - cs.Northf) / s.A
	x2 := x * x
	y2 := y * y
	V1 := 3*x*y2 - x*x2
	V2 := y*y2 - 3*x2*y
	V3 := x * (5*y2*y2 + x2*(-10*y2+x2))
	V4 := y * (5*x2*x2 + y2*(-10*x2+y2))
	x += -q.Ca*V1 - q.Cb*V2 + q.Cc*V3 + q.Cd*V4
	y += q.Cb*V1 - q.Ca*V2 - q.Cd*V3 + q.Cc*V4
	ps := q.p0s + y/q.kRg
	pe := ps + phi0 - q.p0s

	for range 20 {
		V1 = q.A * math.Log(math.Tan(math.Pi/4+0.5*pe))
		tpe := e * math.Sin(pe)
		V2 = 0.5 * e * q.A * math.Log((1+tpe)/(1-tpe))
		t := ps - 2*(math.Atan(math.Exp(V1-V2+q.C))-math.Pi/4)
		pe += t
		if math.Abs(t) < 1e-10 {
			break
		}
	}

	t := e * math.Sin(pe)
	t = 1 - t*t
	Re := (1 - e2) / (t * math.Sqrt(t))
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	tanps := math.Tan(ps)
	t2 := tanps * tanps
	s2 := q.kRg * q.kRg
	d := Re * k0 * q.kRg
	I7 := tanps / (2 * d)
	I8 := tanps * (5 + 3*t2) / (24 * d * s2)
	d = math.Cos(ps) * q.kRg * q.A
	I9 := 1 / d
	d *= s2
	I10 := (1 + 2*t2) / (6 * d)
	I11 := (5 + t2*(28+24*t2)) / (120 * d * s2)
	x2 = x * x
	phi := pe + x2*(-I7+I8*x2)
	lam := x * (I9 + x2*(-I10+x2*I11))

	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
