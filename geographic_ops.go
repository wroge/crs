package crs

type LongitudeRotation struct {
	Lon float64
}

func (l LongitudeRotation) String() string {
	return build("longitude_rotation").addAll("lon", l.Lon).String()
}

func (l LongitudeRotation) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return lon + l.Lon, lat, h, nil
}

func (l LongitudeRotation) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return lon0 - l.Lon, lat0, h0, nil
}

type GeographicOffset struct {
	Lat float64 // arc-seconds
	Lon float64 // arc-seconds
}

func (o GeographicOffset) String() string {
	return build("geographic_offset").addAll(
		"lat", o.Lat,
		"lon", o.Lon,
	).String()
}

func (o GeographicOffset) deg() (dLon, dLat float64) {
	return o.Lon / 3600, o.Lat / 3600
}

func (o GeographicOffset) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	dLon, dLat := o.deg()
	return lon + dLon, lat + dLat, h, nil
}

func (o GeographicOffset) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	dLon, dLat := o.deg()
	return lon0 - dLon, lat0 - dLat, h0, nil
}

type Geographic3dOffset struct {
	Lat float64 // arc-seconds
	Lon float64 // arc-seconds
	H   float64 // metres
}

func (o Geographic3dOffset) String() string {
	return build("geographic_3d_offset").addAll(
		"lat", o.Lat,
		"lon", o.Lon,
		"h", o.H,
	).String()
}

func (o Geographic3dOffset) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return lon + o.Lon/3600, lat + o.Lat/3600, h + o.H, nil
}

func (o Geographic3dOffset) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return lon0 - o.Lon/3600, lat0 - o.Lat/3600, h0 - o.H, nil
}
