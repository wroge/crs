package crs

import (
	"math"
)

type PositionVector struct {
	Tx float64
	Ty float64
	Tz float64
	Rx float64
	Ry float64
	Rz float64
	Ds float64
}

func (pc PositionVector) String() string {
	return build("position_vector").addAll(
		"tx", pc.Tx,
		"ty", pc.Ty,
		"tz", pc.Tz,
		"rx", pc.Rx,
		"ry", pc.Ry,
		"rz", pc.Rz,
		"ds", pc.Ds,
	).String()
}

func (pv PositionVector) FromTarget(source Spheroid, target Spheroid, lon0 float64, lat0 float64, h0 float64) (float64, float64, float64, error) {
	x0, y0, z0 := source.GeographicToGeocentric(lon0, lat0, h0)

	x, y, z := calcHelmert(x0, y0, z0, -pv.Tx, -pv.Ty, -pv.Tz, -pv.Rx, -pv.Ry, -pv.Rz, -pv.Ds)

	lon, lat, h := target.GeocentricToGeographic(x, y, z)

	return lon, lat, h, nil
}

func (pv PositionVector) ToTarget(source Spheroid, target Spheroid, lon float64, lat float64, h float64) (float64, float64, float64, error) {
	x, y, z := source.GeographicToGeocentric(lon, lat, h)

	x0, y0, z0 := calcHelmert(x, y, z, pv.Tx, pv.Ty, pv.Tz, pv.Rx, pv.Ry, pv.Rz, pv.Ds)

	lon0, lat0, h0 := target.GeocentricToGeographic(x0, y0, z0)

	return lon0, lat0, h0, nil
}

func calcHelmert(x, y, z, tx, ty, tz, rx, ry, rz, ds float64) (x1, y1, z1 float64) {
	const (
		asec = math.Pi / 648000
		ppm  = 0.000001
	)

	x1 = (1+ds*ppm)*(x+z*ry*asec-y*rz*asec) + tx
	y1 = (1+ds*ppm)*(y+x*rz*asec-z*rx*asec) + ty
	z1 = (1+ds*ppm)*(z+y*rx*asec-x*ry*asec) + tz

	return x1, y1, z1
}
