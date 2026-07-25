package crs

import "fmt"

type Operation interface {
	ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error)
	FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error)
}

type Transformation struct {
	Target      *Datum
	Accuracy    float64
	BoundingBox BoundingBox
	Operation   Operation
}

func (t Transformation) requireTarget() error {
	if t.Target == nil || t.Target.Name == "" {
		return UnsupportedError{
			Err: fmt.Errorf("no target"),
		}
	}

	return nil
}

func (t Transformation) String() string {
	return build("").addAll(
		"operation", t.Operation,
		"target", t.Target.Name,
		"accuracy", t.Accuracy,
		"bbox", t.BoundingBox,
	).String()
}

func (t Transformation) ToWGS84(source Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return t.toWGS84Visited(source, lon, lat, h, nil, 0)
}

// toWGS84Visited applies this hop toward WGS 84. greenwichLon is lon expressed east
// of Greenwich for area-of-use checks; lon itself may be relative to a local PM.
func (t Transformation) toWGS84Visited(source Spheroid, lon, lat, h float64, visited map[string]bool, greenwichLon float64) (float64, float64, float64, error) {
	if err := t.requireTarget(); err != nil {
		return 0, 0, 0, err
	}

	if !t.BoundingBox.Contains(greenwichLon, lat) {
		return 0, 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("[%f,%f]: %s", greenwichLon, lat, t.BoundingBox),
		}
	}

	lon, lat, h, err := t.Operation.ToTarget(source, t.Target.Spheroid, lon, lat, h)
	if err != nil {
		return 0, 0, 0, err
	}

	return t.Target.toWGS84Visited(lon, lat, h, visited)
}

func (t Transformation) ToDatum(source Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	if err := t.requireTarget(); err != nil {
		return 0, 0, 0, err
	}

	return t.Operation.ToTarget(source, t.Target.Spheroid, lon, lat, h)
}

func (t Transformation) FromDatum(owner Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	if err := t.requireTarget(); err != nil {
		return 0, 0, 0, err
	}

	return t.Operation.FromTarget(t.Target.Spheroid, owner, lon, lat, h)
}

func (t Transformation) FromWGS84(source Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return t.fromWGS84Visited(source, lon0, lat0, h0, nil)
}

func (t Transformation) fromWGS84Visited(source Spheroid, lon0, lat0, h0 float64, visited map[string]bool) (float64, float64, float64, error) {
	if err := t.requireTarget(); err != nil {
		return 0, 0, 0, err
	}

	if !t.BoundingBox.Contains(lon0, lat0) {
		return 0, 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("[%f,%f]: %s", lon0, lat0, t.BoundingBox),
		}
	}

	lon0, lat0, h0, err := t.Target.fromWGS84Visited(lon0, lat0, h0, visited)
	if err != nil {
		return 0, 0, 0, err
	}

	return t.Operation.FromTarget(t.Target.Spheroid, source, lon0, lat0, h0)
}
