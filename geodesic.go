package crs

import "math"

const (
	vincentyEps = 1e-12
	vincentyMax = 100
)

// geodesicInverse is Vincenty's inverse: distance in metres and forward
// azimuth in radians from (lon1, lat1) to (lon2, lat2), degrees.
func geodesicInverse(s Spheroid, lon1, lat1, lon2, lat2 float64) (dist, azi float64, ok bool) {
	if s.A == 0 {
		return 0, 0, false
	}
	if s.F() == 0 {
		return sphereInverse(s.A, lon1, lat1, lon2, lat2)
	}

	a, b, f := s.A, s.B(), s.F()
	u1s, u1c := reducedLat(radian(lat1), f)
	u2s, u2c := reducedLat(radian(lat2), f)
	lonDiff := deltaLon(lon2, lon1)
	lam := lonDiff

	var sinSig, cosSig, sig, sinA, cos2A, cos2Sm float64
	for range vincentyMax {
		sinLam, cosLam := math.Sincos(lam)
		sinSig = math.Hypot(u2c*sinLam, u1c*u2s-u1s*u2c*cosLam)
		if sinSig == 0 {
			return 0, 0, true
		}
		cosSig = u1s*u2s + u1c*u2c*cosLam
		sig = math.Atan2(sinSig, cosSig)
		sinA = u1c * u2c * sinLam / sinSig
		cos2A = 1 - sinA*sinA
		if cos2A == 0 {
			cos2Sm = 0
		} else {
			cos2Sm = cosSig - 2*u1s*u2s/cos2A
		}
		c := f / 16 * cos2A * (4 + f*(4-3*cos2A))
		prev := lam
		lam = lonDiff + (1-c)*f*sinA*(sig+c*sinSig*(cos2Sm+c*cosSig*(-1+2*cos2Sm*cos2Sm)))
		if math.Abs(lam-prev) < vincentyEps {
			u2 := 0.0
			if cos2A != 0 {
				u2 = cos2A * (a*a - b*b) / (b * b)
			}
			aa := 1 + u2/16384*(4096+u2*(-768+u2*(320-175*u2)))
			bb := u2 / 1024 * (256 + u2*(-128+u2*(74-47*u2)))
			dSig := bb * sinSig * (cos2Sm + bb/4*(cosSig*(-1+2*cos2Sm*cos2Sm)-bb/6*cos2Sm*(-3+4*sinSig*sinSig)*(-3+4*cos2Sm*cos2Sm)))
			dist = b * aa * (sig - dSig)
			azi = math.Atan2(u2c*sinLam, u1c*u2s-u1s*u2c*cosLam)
			return dist, azi, true
		}
	}
	return 0, 0, false
}

// geodesicDirect is Vincenty's direct problem. lon/lat are degrees, azi radians,
// dist metres. The returned longitude is in (−180, 180].
func geodesicDirect(s Spheroid, lon1, lat1, azi, dist float64) (lon2, lat2 float64, ok bool) {
	if s.A == 0 {
		return 0, 0, false
	}
	if dist == 0 {
		return lon1, lat1, true
	}
	if s.F() == 0 {
		return sphereDirect(s.A, lon1, lat1, azi, dist)
	}

	a, b, f := s.A, s.B(), s.F()
	u1s, u1c := reducedLat(radian(lat1), f)
	sinAzi, cosAzi := math.Sincos(azi)
	sig1 := math.Atan2(u1s, u1c*cosAzi)
	sinA := u1c * sinAzi
	cos2A := 1 - sinA*sinA
	u2 := 0.0
	if cos2A != 0 {
		u2 = cos2A * (a*a - b*b) / (b * b)
	}
	aa := 1 + u2/16384*(4096+u2*(-768+u2*(320-175*u2)))
	bb := u2 / 1024 * (256 + u2*(-128+u2*(74-47*u2)))
	if aa == 0 || b == 0 {
		return 0, 0, false
	}

	sig := dist / (b * aa)
	var cos2Sm, sinSig, cosSig float64
	converged := false
	for range vincentyMax {
		cos2Sm = math.Cos(2*sig1 + sig)
		sinSig, cosSig = math.Sincos(sig)
		dSig := bb * sinSig * (cos2Sm + bb/4*(cosSig*(-1+2*cos2Sm*cos2Sm)-bb/6*cos2Sm*(-3+4*sinSig*sinSig)*(-3+4*cos2Sm*cos2Sm)))
		next := dist/(b*aa) + dSig
		if math.Abs(next-sig) < vincentyEps {
			sig = next
			converged = true
			break
		}
		sig = next
	}
	if !converged {
		return 0, 0, false
	}
	sinSig, cosSig = math.Sincos(sig)

	lat2 = degree(math.Atan2(u1s*cosSig+u1c*sinSig*cosAzi, (1-f)*math.Hypot(sinA, u1s*sinSig-u1c*cosSig*cosAzi)))
	lam := math.Atan2(sinSig*sinAzi, u1c*cosSig-u1s*sinSig*cosAzi)
	c := f / 16 * cos2A * (4 + f*(4-3*cos2A))
	lonDiff := lam - (1-c)*f*sinA*(sig+c*sinSig*(cos2Sm+c*cosSig*(-1+2*cos2Sm*cos2Sm)))
	lon2 = math.Remainder(lon1+degree(lonDiff), 360)
	return lon2, lat2, true
}

func reducedLat(phi, f float64) (sinU, cosU float64) {
	if phi >= math.Pi/2-1e-14 {
		return 1, 0
	}
	if phi <= -math.Pi/2+1e-14 {
		return -1, 0
	}
	sinPhi, cosPhi := math.Sincos(phi)
	y := (1 - f) * sinPhi
	h := math.Hypot(y, cosPhi)
	if h == 0 {
		return 0, 1
	}
	return y / h, cosPhi / h
}

func sphereInverse(a, lon1, lat1, lon2, lat2 float64) (dist, azi float64, ok bool) {
	if a == 0 {
		return 0, 0, false
	}
	dLam := deltaLon(lon2, lon1)
	sin1, cos1 := math.Sincos(radian(lat1))
	sin2, cos2 := math.Sincos(radian(lat2))
	sinLam, cosLam := math.Sincos(dLam)
	cosC := sin1*sin2 + cos1*cos2*cosLam
	c := math.Acos(clamp(cosC, -1, 1))
	if c == 0 {
		return 0, 0, true
	}
	return a * c, math.Atan2(cos2*sinLam, cos1*sin2-sin1*cos2*cosLam), true
}

func sphereDirect(a, lon1, lat1, azi, dist float64) (lon2, lat2 float64, ok bool) {
	if a == 0 {
		return 0, 0, false
	}
	c := dist / a
	sinC, cosC := math.Sincos(c)
	sin1, cos1 := math.Sincos(radian(lat1))
	sinAzi, cosAzi := math.Sincos(azi)
	lat2 = degree(math.Asin(clamp(sin1*cosC+cos1*sinC*cosAzi, -1, 1)))
	lam := math.Atan2(sinAzi*sinC, cos1*cosC-sin1*sinC*cosAzi)
	lon2 = math.Remainder(lon1+degree(lam), 360)
	return lon2, lat2, true
}
