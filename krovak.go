package crs

import (
	"math"
)

type Krovak struct {
	Lonf   float64
	Latf   float64
	Alpha  float64
	Scale  float64
	Eastf  float64
	Northf float64
}

func (cs Krovak) String() string {
	return build("krovak").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"alpha", cs.Alpha,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (k Krovak) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phic := radian(k.Latf)
	lambda0 := radian(k.Lonf)
	phip := radian(78.5)
	alphac := radian(k.Alpha)

	A := s.A * math.Sqrt(1-s.E2()) / (1 - s.E2()*sin2(phic))
	B := math.Sqrt(1 + (s.E2() * intPow(math.Cos(phic), 4) / (1 - s.E2())))
	gamma0 := math.Asin(math.Sin(phic) / B)
	t0 := math.Tan(math.Pi/4+gamma0/2) * math.Pow((1+s.E()*math.Sin(phic))/(1-s.E()*math.Sin(phic)), s.E()*B/2) / math.Pow(math.Tan(math.Pi/4+phic/2), B)
	n := math.Sin(phip)
	r0 := k.Scale * A / math.Tan(phip)

	phi := radian(lat)
	lambda := radian(lon)

	U := 2 * (math.Atan(t0*math.Pow(math.Tan(phi/2+math.Pi/4), B)/math.Pow((1+s.E()*math.Sin(phi))/(1-s.E()*math.Sin(phi)), s.E()*B/2)) - math.Pi/4)
	V := B * (lambda0 - lambda)
	T := math.Asin(math.Cos(alphac)*math.Sin(U) + math.Sin(alphac)*math.Cos(U)*math.Cos(V))
	D := math.Asin(math.Cos(U) * math.Sin(V) / math.Cos(T))
	theta := n * D
	r := r0 * math.Pow(math.Tan(math.Pi/4+phip/2), n) / math.Pow(math.Tan(T/2+math.Pi/4), n)
	Xp := r * math.Cos(theta)
	Yp := r * math.Sin(theta)

	return -(Yp + k.Eastf), -(Xp + k.Northf), h
}

func (k Krovak) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phic := radian(k.Latf)
	lambda0 := radian(k.Lonf)
	phip := radian(78.5)
	alphac := radian(k.Alpha)

	A := s.A * math.Sqrt(1-s.E2()) / (1 - s.E2()*sin2(phic))
	B := math.Sqrt(1 + (s.E2() * intPow(math.Cos(phic), 4) / (1 - s.E2())))
	gamma0 := math.Asin(math.Sin(phic) / B)
	t0 := math.Tan(math.Pi/4+gamma0/2) * math.Pow((1+s.E()*math.Sin(phic))/(1-s.E()*math.Sin(phic)), s.E()*B/2) / math.Pow(math.Tan(math.Pi/4+phic/2), B)
	n := math.Sin(phip)
	r0 := k.Scale * A / math.Tan(phip)

	Xpi := (-north) - k.Northf
	Ypi := (-east) - k.Eastf
	ri := math.Sqrt(intPow(Xpi, 2) + intPow(Ypi, 2))
	thetai := math.Atan2(Ypi, Xpi)
	di := thetai / math.Sin(phip)
	ti := 2 * (math.Atan(math.Pow(r0/ri, 1/n)*math.Tan(math.Pi/4+phip/2)) - math.Pi/4)
	ui := math.Asin(math.Cos(alphac)*math.Sin(ti) - math.Sin(alphac)*math.Cos(ti)*math.Cos(di))
	vi := math.Asin(math.Cos(ti) * math.Sin(di) / math.Cos(ui))

	phi := ui

	for range 12 {
		phiNext := 2 * (math.Atan(math.Pow(t0, -1/B)*math.Pow(math.Tan(ui/2+math.Pi/4), 1/B)*math.Pow((1+s.E()*math.Sin(phi))/(1-s.E()*math.Sin(phi)), s.E()/2)) - math.Pi/4)
		if math.Abs(phiNext-phi) < 1e-12 {
			phi = phiNext
			break
		}
		phi = phiNext
	}

	lambda := lambda0 - vi/B

	return degree(lambda), degree(phi), h
}
