package crs

import (
	"math"
)

const (
	krovakModX0 = 1089000.0
	krovakModY0 = 654000.0
)

const (
	krovakModC1  = 2.946529277e-02
	krovakModC2  = 2.515965696e-02
	krovakModC3  = 1.193845912e-07
	krovakModC4  = -4.668270147e-07
	krovakModC5  = 9.233980362e-12
	krovakModC6  = 1.523735715e-12
	krovakModC7  = 1.696780024e-18
	krovakModC8  = 4.408314235e-18
	krovakModC9  = -8.331083518e-24
	krovakModC10 = -3.689471323e-24
)

type KrovakModified struct {
	Lonf, Latf, Alpha, Scale, Eastf, Northf float64
}

func (cs KrovakModified) String() string {
	return build("krovak_modified").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"alpha", cs.Alpha,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs KrovakModified) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	xp, yp := krovakXY(s, lon, lat, cs.Lonf, cs.Latf, cs.Alpha, cs.Scale)
	xp, yp = krovakModifiedApply(xp, yp)

	return -(yp + cs.Eastf), -(xp + cs.Northf), h
}

func (cs KrovakModified) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	xp, yp := krovakModifiedInvert((-north)-cs.Northf, (-east)-cs.Eastf)

	return krovakLonLat(s, xp, yp, cs.Lonf, cs.Latf, cs.Alpha, cs.Scale, h)
}

func krovakXY(s Spheroid, lon, lat, lonf, latf, alpha, scale float64) (xp, yp float64) {
	phic := radian(latf)
	lambda0 := radian(lonf)
	phip := radian(78.5)
	alphac := radian(alpha)

	A := s.A * math.Sqrt(1-s.E2()) / (1 - s.E2()*sin2(phic))
	B := math.Sqrt(1 + (s.E2() * intPow(math.Cos(phic), 4) / (1 - s.E2())))
	gamma0 := math.Asin(math.Sin(phic) / B)
	t0 := math.Tan(math.Pi/4+gamma0/2) * math.Pow((1+s.E()*math.Sin(phic))/(1-s.E()*math.Sin(phic)), s.E()*B/2) / math.Pow(math.Tan(math.Pi/4+phic/2), B)
	n := math.Sin(phip)
	r0 := scale * A / math.Tan(phip)

	phi := radian(lat)
	lambda := radian(lon)

	U := 2 * (math.Atan(t0*math.Pow(math.Tan(phi/2+math.Pi/4), B)/math.Pow((1+s.E()*math.Sin(phi))/(1-s.E()*math.Sin(phi)), s.E()*B/2)) - math.Pi/4)
	V := B * (lambda0 - lambda)
	T := math.Asin(math.Cos(alphac)*math.Sin(U) + math.Sin(alphac)*math.Cos(U)*math.Cos(V))
	D := math.Asin(math.Cos(U) * math.Sin(V) / math.Cos(T))
	theta := n * D
	r := r0 * math.Pow(math.Tan(math.Pi/4+phip/2), n) / math.Pow(math.Tan(T/2+math.Pi/4), n)

	return r * math.Cos(theta), r * math.Sin(theta)
}

func krovakLonLat(s Spheroid, xp, yp, lonf, latf, alpha, scale, h float64) (float64, float64, float64) {
	phic := radian(latf)
	lambda0 := radian(lonf)
	phip := radian(78.5)
	alphac := radian(alpha)

	A := s.A * math.Sqrt(1-s.E2()) / (1 - s.E2()*sin2(phic))
	B := math.Sqrt(1 + (s.E2() * intPow(math.Cos(phic), 4) / (1 - s.E2())))
	gamma0 := math.Asin(math.Sin(phic) / B)
	t0 := math.Tan(math.Pi/4+gamma0/2) * math.Pow((1+s.E()*math.Sin(phic))/(1-s.E()*math.Sin(phic)), s.E()*B/2) / math.Pow(math.Tan(math.Pi/4+phic/2), B)
	n := math.Sin(phip)
	r0 := scale * A / math.Tan(phip)

	ri := math.Sqrt(intPow(xp, 2) + intPow(yp, 2))
	thetai := math.Atan2(yp, xp)
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

	return degree(lambda0 - vi/B), degree(phi), h
}

func krovakModifiedDXDY(xr, yr float64) (dX, dY float64) {
	xr2 := xr * xr
	yr2 := yr * yr
	xr4 := xr2 * xr2
	yr4 := yr2 * yr2
	dX = krovakModC1 + krovakModC3*xr - krovakModC4*yr - 2*krovakModC6*xr*yr + krovakModC5*(xr2-yr2) +
		krovakModC7*xr*(xr2-3*yr2) - krovakModC8*yr*(3*xr2-yr2) +
		4*krovakModC9*xr*yr*(xr2-yr2) + krovakModC10*(xr4+yr4-6*xr2*yr2)
	dY = krovakModC2 + krovakModC3*yr + krovakModC4*xr + 2*krovakModC5*xr*yr + krovakModC6*(xr2-yr2) +
		krovakModC8*xr*(xr2-3*yr2) + krovakModC7*yr*(3*xr2-yr2) -
		4*krovakModC10*xr*yr*(xr2-yr2) + krovakModC9*(xr4+yr4-6*xr2*yr2)

	return dX, dY
}

func krovakModifiedApply(xp, yp float64) (float64, float64) {
	dX, dY := krovakModifiedDXDY(xp-krovakModX0, yp-krovakModY0)

	return xp - dX, yp - dY
}

func krovakModifiedInvert(xp, yp float64) (float64, float64) {
	u, v := xp, yp

	for range 10 {
		dX, dY := krovakModifiedDXDY(u-krovakModX0, v-krovakModY0)
		nu, nv := xp+dX, yp+dY
		if math.Abs(nu-u) < 1e-12 && math.Abs(nv-v) < 1e-12 {
			return nu, nv
		}
		u, v = nu, nv
	}

	return u, v
}
