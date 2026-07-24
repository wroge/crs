package crs

import "math"

type CoordinateFrameFullMatrix struct {
	Tx, Ty, Tz, Rx, Ry, Rz, Ds float64
}

func (m CoordinateFrameFullMatrix) String() string {
	return build("coordinate_frame_full_matrix").addAll(
		"tx", m.Tx,
		"ty", m.Ty,
		"tz", m.Tz,
		"rx", m.Rx,
		"ry", m.Ry,
		"rz", m.Rz,
		"ds", m.Ds,
	).String()
}

func (m CoordinateFrameFullMatrix) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	x, y, z := source.GeographicToGeocentric(lon, lat, h)
	x0, y0, z0 := calcHelmertExact(x, y, z, m.Tx, m.Ty, m.Tz, m.Rx, m.Ry, m.Rz, m.Ds)
	lon0, lat0, h0 := target.GeocentricToGeographic(x0, y0, z0)

	return lon0, lat0, h0, nil
}

func (m CoordinateFrameFullMatrix) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	x0, y0, z0 := source.GeographicToGeocentric(lon0, lat0, h0)
	x, y, z := calcHelmertExactInverse(x0, y0, z0, m.Tx, m.Ty, m.Tz, m.Rx, m.Ry, m.Rz, m.Ds)
	lon, lat, h := target.GeocentricToGeographic(x, y, z)

	return lon, lat, h, nil
}

const (
	asec = math.Pi / 648000
	ppm  = 0.000001
)

func calcHelmertExact(x, y, z, tx, ty, tz, rx, ry, rz, ds float64) (x1, y1, z1 float64) {
	r00, r01, r02, r10, r11, r12, r20, r21, r22 := rotationMatrixExact(rx*asec, ry*asec, rz*asec)
	s := 1 + ds*ppm
	x1 = s*(r00*x+r01*y+r02*z) + tx
	y1 = s*(r10*x+r11*y+r12*z) + ty
	z1 = s*(r20*x+r21*y+r22*z) + tz

	return x1, y1, z1
}

func calcHelmertExactInverse(x1, y1, z1, tx, ty, tz, rx, ry, rz, ds float64) (x, y, z float64) {
	r00, r01, r02, r10, r11, r12, r20, r21, r22 := rotationMatrixExact(rx*asec, ry*asec, rz*asec)

	s := 1 + ds*ppm
	dx, dy, dz := x1-tx, y1-ty, z1-tz
	x = (r00*dx + r10*dy + r20*dz) / s
	y = (r01*dx + r11*dy + r21*dz) / s
	z = (r02*dx + r12*dy + r22*dz) / s

	return x, y, z
}

func rotationMatrixExact(rx, ry, rz float64) (r00, r01, r02, r10, r11, r12, r20, r21, r22 float64) {
	sx, cx := math.Sincos(rx)
	sy, cy := math.Sincos(ry)
	sz, cz := math.Sincos(rz)

	r00 = cy * cz
	r01 = -cx*sz + sx*sy*cz
	r02 = sx*sz + cx*sy*cz
	r10 = cy * sz
	r11 = cx*cz + sx*sy*sz
	r12 = -sx*cz + cx*sy*sz
	r20 = -sy
	r21 = sx * cy
	r22 = cx * cy
	return
}
