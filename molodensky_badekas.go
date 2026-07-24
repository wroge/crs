package crs

type MolodenskyBadekas struct {
	Tx, Ty, Tz, Rx, Ry, Rz, Ds float64
	Px, Py, Pz                 float64
}

func (m MolodenskyBadekas) String() string {
	return build("molodensky_badekas").addAll(
		"tx", m.Tx,
		"ty", m.Ty,
		"tz", m.Tz,
		"rx", m.Rx,
		"ry", m.Ry,
		"rz", m.Rz,
		"ds", m.Ds,
		"px", m.Px,
		"py", m.Py,
		"pz", m.Pz,
	).String()
}

func (m MolodenskyBadekas) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	x, y, z := source.GeographicToGeocentric(lon, lat, h)
	x0, y0, z0 := calcMolodenskyBadekas(x, y, z, m.Tx, m.Ty, m.Tz, m.Rx, m.Ry, m.Rz, m.Ds, m.Px, m.Py, m.Pz)
	lon0, lat0, h0 := target.GeocentricToGeographic(x0, y0, z0)
	return lon0, lat0, h0, nil
}

func (m MolodenskyBadekas) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	x0, y0, z0 := source.GeographicToGeocentric(lon0, lat0, h0)
	x, y, z := calcMolodenskyBadekasInverse(x0, y0, z0, m.Tx, m.Ty, m.Tz, m.Rx, m.Ry, m.Rz, m.Ds, m.Px, m.Py, m.Pz)
	lon, lat, h := target.GeocentricToGeographic(x, y, z)

	return lon, lat, h, nil
}

func calcMolodenskyBadekas(x, y, z, tx, ty, tz, rx, ry, rz, ds, px, py, pz float64) (x1, y1, z1 float64) {
	x1, y1, z1 = calcHelmert(x-px, y-py, z-pz, tx, ty, tz, rx, ry, rz, ds)

	return x1 + px, y1 + py, z1 + pz
}

func calcMolodenskyBadekasInverse(x1, y1, z1, tx, ty, tz, rx, ry, rz, ds, px, py, pz float64) (x, y, z float64) {
	x, y, z = calcHelmert(x1-tx-px, y1-ty-py, z1-tz-pz, 0, 0, 0, -rx, -ry, -rz, -ds)
	return x + px, y + py, z + pz
}
