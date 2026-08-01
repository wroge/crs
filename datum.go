package crs

import (
	"embed"
	"errors"
	"fmt"
	"io"
	"slices"
	"strconv"
	"strings"
	"sync"
)

var (
	WGS84 = Datum{
		Name: "wgs84",
		Spheroid: Spheroid{
			Name: "wgs84",
			A:    6378137,
			Fi:   298.257223563,
		},
	}
)

//go:embed datum/*.txt
var datumDir embed.FS

var datumStore sync.Map

func RegisterDatum(d Datum) {
	datumStore.Store(d.Name, d)
}

func loadDatum(name string) (d Datum, err error) {
	name = strings.ToLower(name)

	v, ok := datumStore.Load(name)
	if ok {
		return v.(Datum), nil
	}

	switch name {
	case "wgs84":
		return WGS84, nil
	}

	defer func() {
		if err == nil {
			datumStore.Store(name, d)
		}
	}()

	file, err := datumDir.Open(fmt.Sprintf("datum/%s.txt", name))
	if err != nil {
		return d, fmt.Errorf("datum not found: %s", name)
	}

	data, err := io.ReadAll(file)
	if err != nil {
		return d, err
	}

	d, err = parseDatum(string(data))
	if err != nil {
		return d, err
	}

	d.Name = name

	return d, nil
}

func stripComment(line string) string {
	if i := strings.Index(line, "#"); i >= 0 {
		line = line[:i]
	}

	return strings.TrimSpace(line)
}

func asParts(txt string) parts {
	txt = flattenDSL(txt)
	var pp []part

	for f := range strings.FieldsSeq(txt) {
		key, value, ok := strings.Cut(f, "=")
		if ok {
			pp = append(pp, part{
				key:   strings.Trim(key, "-+"),
				value: value,
			})

			continue
		}

		key, value, ok = strings.Cut(f, ":")
		if ok {
			pp = append(pp, part{
				key:   strings.Trim(key, "-+"),
				value: value,
			})

			continue
		}

		pp = append(pp, part{
			key: strings.Trim(key, "-+"),
		})
	}

	return pp
}

type parts []part

func (pp parts) getPart(find ...string) part {
	for _, part := range pp {
		if slices.Contains(find, part.key) {
			return part
		}
	}

	return part{}
}

func (pp parts) asTransformation() (Transformation, error) {
	var t Transformation

	t.Accuracy, _ = pp.getPart("accuracy").asFloat()

	bbox, ok := pp.getPart("bbox").asBoundingBox()
	if ok {
		t.BoundingBox = bbox
	} else {
		t.BoundingBox = World
	}

	targetName, ok := pp.getPart("target").asString()
	if !ok {
		targetName = "wgs84"
	}

	target, err := loadDatum(targetName)
	if err != nil {
		if _, serr := loadSpheroid(targetName); serr == nil {
			return Transformation{}, fmt.Errorf("target %q is a spheroid, not a datum", targetName)
		}
		return Transformation{}, err
	}

	if target.Name == "" {
		return Transformation{}, UnsupportedError{
			Err: errors.New("no target"),
		}
	}
	t.Target = &target

	op, err := pp.asOperation()
	if err != nil {
		if errors.As(err, &UnsupportedError{}) {
			return t, nil
		}

		return Transformation{}, err
	}

	t.Operation = op

	return t, nil
}

type part struct {
	key   string
	value string
	err   error
}

func (p part) asString() (string, bool) {
	if p.err != nil || p.value == "" {
		return "", false
	}

	return p.value, true
}

func (p part) asFloat() (float64, bool) {
	if p.err != nil || p.value == "" {
		return 0, false
	}

	v, err := strconv.ParseFloat(p.value, 64)
	if err != nil {
		return 0, false
	}

	return v, true
}

func (pp parts) float(key string) (float64, bool) {
	return pp.getPart(key).asFloat()
}

func (pp parts) floatOrZero(key string) float64 {
	v, _ := pp.float(key)
	return v
}

func (pp parts) hasInverse() bool {
	for _, p := range pp {
		if p.key == "inverse" {
			return true
		}
	}
	return false
}

func (pp parts) asOperation() (Operation, error) {
	opName, ok := pp.getPart("operation").asString()
	if !ok {
		return nil, nil
	}

	var op Operation

	switch opName {
	case "position_vector":
		op = PositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: pp.floatOrZero("rx"), Ry: pp.floatOrZero("ry"), Rz: pp.floatOrZero("rz"),
			Ds: pp.floatOrZero("ds"),
		}
	case "coordinate_frame":
		op = PositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: -pp.floatOrZero("rx"), Ry: -pp.floatOrZero("ry"), Rz: -pp.floatOrZero("rz"),
			Ds: pp.floatOrZero("ds"),
		}
	case "horizontal_grid":
		grid, ok := pp.getPart("grid").asString()
		if !ok {
			return nil, fmt.Errorf("horizontal_grid missing grid=")
		}
		op = HorizontalGrid(grid)
	case "vertical_grid":
		grid, ok := pp.getPart("grid").asString()
		if !ok {
			return nil, fmt.Errorf("vertical_grid missing grid=")
		}
		op = VerticalGrid(grid)
	case "vertical_offset":
		op = VerticalOffset{Dh: pp.floatOrZero("dh")}
	case "vertical_offset_and_slope":
		op = VerticalOffsetAndSlope{
			Lat0:     pp.floatOrZero("lat0"),
			Lon0:     pp.floatOrZero("lon0"),
			Dh:       pp.floatOrZero("dh"),
			SlopeLat: pp.floatOrZero("slope_lat"),
			SlopeLon: pp.floatOrZero("slope_lon"),
		}
	case "geocentric_translations":
		op = PositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
		}
	case "identity":
		op = Identity{}
	case "longitude_rotation":
		op = LongitudeRotation{Lon: pp.floatOrZero("lon")}
	case "geographic_offset":
		op = GeographicOffset{
			Lat: pp.floatOrZero("lat"),
			Lon: pp.floatOrZero("lon"),
		}
	case "geographic_3d_offset":
		op = Geographic3dOffset{
			Lat: pp.floatOrZero("lat"),
			Lon: pp.floatOrZero("lon"),
			H:   pp.floatOrZero("h"),
		}
	case "molodensky_badekas":
		op = MolodenskyBadekas{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: -pp.floatOrZero("rx"), Ry: -pp.floatOrZero("ry"), Rz: -pp.floatOrZero("rz"),
			Ds: pp.floatOrZero("ds"),
			Px: pp.floatOrZero("px"), Py: pp.floatOrZero("py"), Pz: pp.floatOrZero("pz"),
		}
	case "molodensky_badekas_pv", "molodensky_badekas_pv_geocentric":
		op = MolodenskyBadekas{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: pp.floatOrZero("rx"), Ry: pp.floatOrZero("ry"), Rz: pp.floatOrZero("rz"),
			Ds: pp.floatOrZero("ds"),
			Px: pp.floatOrZero("px"), Py: pp.floatOrZero("py"), Pz: pp.floatOrZero("pz"),
		}
	case "time_specific_position_vector":
		op = TimeSpecificPositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: pp.floatOrZero("rx"), Ry: pp.floatOrZero("ry"), Rz: pp.floatOrZero("rz"),
			Ds:                  pp.floatOrZero("ds"),
			TransformationEpoch: pp.floatOrZero("epoch"),
		}
	case "time_specific_coordinate_frame":
		op = TimeSpecificPositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: -pp.floatOrZero("rx"), Ry: -pp.floatOrZero("ry"), Rz: -pp.floatOrZero("rz"),
			Ds:                  pp.floatOrZero("ds"),
			TransformationEpoch: pp.floatOrZero("epoch"),
		}
	case "time_dependent_position_vector":
		op = TimeDependentPositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: pp.floatOrZero("rx"), Ry: pp.floatOrZero("ry"), Rz: pp.floatOrZero("rz"),
			Ds:             pp.floatOrZero("ds"),
			TxRate:         pp.floatOrZero("dtx"),
			TyRate:         pp.floatOrZero("dty"),
			TzRate:         pp.floatOrZero("dtz"),
			RxRate:         pp.floatOrZero("drx"),
			RyRate:         pp.floatOrZero("dry"),
			RzRate:         pp.floatOrZero("drz"),
			DsRate:         pp.floatOrZero("dds"),
			ReferenceEpoch: pp.floatOrZero("epoch"),
		}
	case "time_dependent_coordinate_frame":
		op = TimeDependentPositionVector{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: -pp.floatOrZero("rx"), Ry: -pp.floatOrZero("ry"), Rz: -pp.floatOrZero("rz"),
			Ds:             pp.floatOrZero("ds"),
			TxRate:         pp.floatOrZero("dtx"),
			TyRate:         pp.floatOrZero("dty"),
			TzRate:         pp.floatOrZero("dtz"),
			RxRate:         -pp.floatOrZero("drx"),
			RyRate:         -pp.floatOrZero("dry"),
			RzRate:         -pp.floatOrZero("drz"),
			DsRate:         pp.floatOrZero("dds"),
			ReferenceEpoch: pp.floatOrZero("epoch"),
		}
	case "coordinate_frame_full_matrix":
		op = CoordinateFrameFullMatrix{
			Tx: pp.floatOrZero("tx"), Ty: pp.floatOrZero("ty"), Tz: pp.floatOrZero("tz"),
			Rx: -pp.floatOrZero("rx"), Ry: -pp.floatOrZero("ry"), Rz: -pp.floatOrZero("rz"),
			Ds: pp.floatOrZero("ds"),
		}
	case "velocity_grid":
		grid, ok := pp.getPart("grid").asString()
		if !ok {
			return nil, fmt.Errorf("velocity_grid missing grid=")
		}
		op = VelocityGrid{
			Grid:  grid,
			Dt:    pp.floatOrZero("dt"),
			Epoch: pp.floatOrZero("epoch"),
		}
	default:
		return nil, UnsupportedError{
			Err: fmt.Errorf("invalid operation: %s", opName),
		}
	}

	if pp.hasInverse() {
		return Inverse{
			Operation: op,
		}, nil
	}

	return op, nil
}

func (p part) asBoundingBox() (BoundingBox, bool) {
	var (
		b   BoundingBox
		err error
	)

	if p.err != nil || p.value == "" {
		return b, false
	}

	for i, v := range strings.Split(p.value, ",") {
		if p.value == "" {
			continue
		}

		switch i {
		case 0:
			b.MinLon, err = strconv.ParseFloat(v, 64)
			if err != nil {
				return b, false
			}
		case 1:
			b.MinLat, err = strconv.ParseFloat(v, 64)
			if err != nil {
				return b, false
			}
		case 2:
			b.MaxLon, err = strconv.ParseFloat(v, 64)
			if err != nil {
				return b, false
			}
		case 3:
			b.MaxLat, err = strconv.ParseFloat(v, 64)
			if err != nil {
				return b, false
			}
		}
	}

	return b, true
}

type Datum struct {
	Name            string
	Spheroid        Spheroid
	Transformations []Transformation

	// coordinateEpoch is set by AtEpoch / TransformAt. When hasCoordinateEpoch
	// is false (plain Transform), time-specific Helmerts are excluded from paths.
	coordinateEpoch    float64
	hasCoordinateEpoch bool
}

func (d Datum) String() string {
	return build("").add("spheroid", d.Spheroid.Name).add("", d.Transformations).String()
}

func (d Datum) Intersects(bbox ...BoundingBox) (Datum, error) {
	if len(bbox) == 0 {
		return d, nil
	}

	if len(d.Transformations) == 0 {
		return d, nil
	}

	// Clone so DeleteFunc cannot mutate the cached LoadDatum slice backing array.
	t := slices.DeleteFunc(slices.Clone(d.Transformations), func(t Transformation) bool {
		for _, b := range bbox {
			if !t.BoundingBox.Intersects(b) {
				return true
			}
		}

		return false
	})

	if len(t) == 0 {
		// No published hub ops cover this area (e.g. Guam1963 on Yap). Keep the
		// spheroid so same-datum projection still works; cross-datum will fail later.
		return Datum{
			Name:               d.Name,
			Spheroid:           d.Spheroid,
			coordinateEpoch:    d.coordinateEpoch,
			hasCoordinateEpoch: d.hasCoordinateEpoch,
		}, nil
	}

	return Datum{
		Name:               d.Name,
		Spheroid:           d.Spheroid,
		Transformations:    t,
		coordinateEpoch:    d.coordinateEpoch,
		hasCoordinateEpoch: d.hasCoordinateEpoch,
	}, nil
}

// AtEpoch returns a copy of d with time-dependent operations evaluated at epoch
// (decimal year). TimeDependentPositionVector becomes PositionVector; VelocityGrid
// with a non-zero Epoch gets Dt = epoch - Epoch. Time-specific Helmerts that do
// not match epoch are dropped from the in-memory list (pathfinding also gates them).
func (d Datum) AtEpoch(epoch float64) Datum {
	out := Datum{
		Name:               d.Name,
		Spheroid:           d.Spheroid,
		coordinateEpoch:    epoch,
		hasCoordinateEpoch: true,
	}

	if len(d.Transformations) == 0 {
		return out
	}

	ts := make([]Transformation, 0, len(d.Transformations))
	for _, t := range d.Transformations {
		op := operationAtEpoch(t.Operation, epoch)
		if !timeSpecificAllowed(op, true, epoch) {
			continue
		}

		t.Operation = op
		ts = append(ts, t)
	}

	out.Transformations = ts

	return out
}

func operationAtEpoch(op Operation, epoch float64) Operation {
	switch o := op.(type) {
	case TimeDependentPositionVector:
		return o.at(epoch)
	case Inverse:
		return Inverse{Operation: operationAtEpoch(o.Operation, epoch)}
	case VelocityGrid:
		if o.Epoch != 0 {
			o.Dt = epoch - o.Epoch
		}

		return o
	default:
		return op
	}
}

// primeMeridianLongitude is degrees east of Greenwich for this datum's prime
// meridian, taken from a longitude_rotation hop (0 when the datum is Greenwich-based).
// EPSG areas of use are Greenwich-oriented; use lon+primeMeridianLongitude() for Contains.
func (d Datum) primeMeridianLongitude() float64 {
	for _, t := range d.Transformations {
		if lr, ok := t.Operation.(LongitudeRotation); ok {
			return lr.Lon
		}
	}

	return 0
}

func (d Datum) ToWGS84(lon, lat, h float64) (float64, float64, float64, error) {
	return d.toWGS84Visited(lon, lat, h, nil)
}

func (d Datum) toWGS84Visited(lon, lat, h float64, visited map[string]bool) (float64, float64, float64, error) {
	if len(d.Transformations) == 0 {
		return lon, lat, h, nil
	}

	if visited == nil {
		visited = make(map[string]bool)
	}

	if visited[d.Name] {
		return 0, 0, 0, UnsupportedError{
			Err: errors.New("datum transformation cycle"),
		}
	}

	visited[d.Name] = true

	defer delete(visited, d.Name)

	glon := lon + d.primeMeridianLongitude()
	var lastErr error
	for _, t := range d.Transformations {
		if !timeSpecificAllowed(t.Operation, d.hasCoordinateEpoch, d.coordinateEpoch) {
			continue
		}

		lon0, lat0, h0, err := t.toWGS84Visited(d.Spheroid, lon, lat, h, visited, glon)
		if err != nil {
			lastErr = err
			continue
		}

		return lon0, lat0, h0, nil
	}

	if lastErr != nil {
		return 0, 0, 0, lastErr
	}

	return 0, 0, 0, UnsupportedError{
		Err: errors.New("no valid transformation found"),
	}
}

func (d Datum) FromWGS84(lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return d.fromWGS84Visited(lon0, lat0, h0, nil)
}

func (d Datum) fromWGS84Visited(lon0, lat0, h0 float64, visited map[string]bool) (float64, float64, float64, error) {
	if len(d.Transformations) == 0 {
		return lon0, lat0, h0, nil
	}

	if visited == nil {
		visited = make(map[string]bool)
	}

	if visited[d.Name] {
		return 0, 0, 0, UnsupportedError{
			Err: errors.New("datum transformation cycle"),
		}
	}

	visited[d.Name] = true

	defer delete(visited, d.Name)

	var lastErr error
	for _, t := range d.Transformations {
		if !timeSpecificAllowed(t.Operation, d.hasCoordinateEpoch, d.coordinateEpoch) {
			continue
		}

		lon, lat, h, err := t.fromWGS84Visited(d.Spheroid, lon0, lat0, h0, visited)
		if err != nil {
			lastErr = err
			continue
		}

		return lon, lat, h, nil
	}

	if lastErr != nil {
		return 0, 0, 0, lastErr
	}

	return 0, 0, 0, UnsupportedError{
		Err: errors.New("no valid transformation found"),
	}
}

func (d Datum) TransformTo(target Datum, lon, lat, h float64) (float64, float64, float64, error) {
	if d.Name == target.Name {
		return lon, lat, h, nil
	}

	hasEpoch := d.hasCoordinateEpoch || target.hasCoordinateEpoch
	epoch := d.coordinateEpoch
	if target.hasCoordinateEpoch {
		epoch = target.coordinateEpoch
	}

	excluded := make(map[edgeKey]bool)
	var lastErr error

	for {
		path, err := findBestPath(d, target, lon, lat, excluded, hasEpoch, epoch)
		if err != nil {
			if lastErr != nil {
				return 0, 0, 0, lastErr
			}
			return 0, 0, 0, err
		}

		outLon, outLat, outH, failed, err := applyPath(path, lon, lat, h, excluded, hasEpoch, epoch)
		if err == nil {
			return outLon, outLat, outH, nil
		}

		lastErr = err
		if len(failed) == 0 {
			return 0, 0, 0, lastErr
		}

		grew := false
		for _, k := range failed {
			if !excluded[k] {
				excluded[k] = true
				grew = true
			}
		}

		if !grew {
			return 0, 0, 0, lastErr
		}
	}
}
