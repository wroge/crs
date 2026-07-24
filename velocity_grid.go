package crs

import (
	"errors"
	"fmt"
	"io"
	"math"
	"path"
	"strings"
	"sync"
)

// VelocityGrid is PROJ +proj=deformation with a fixed timespan dt (years).
// Grid bands are east/north/up velocities in mm/year (GTG TYPE=VELOCITY).
// When Epoch != 0, Datum.AtEpoch sets Dt = epoch - Epoch.
type VelocityGrid struct {
	Grid  string
	Dt    float64
	Epoch float64
}

func (v VelocityGrid) String() string {
	return build("velocity_grid").addAll(
		"grid", v.Grid,
		"dt", v.Dt,
		"epoch", v.Epoch,
	).String()
}

// ToTarget implements [Operation]: apply +dt · V in geocentric space.
func (v VelocityGrid) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return v.apply(source, target, lon, lat, h, v.Dt)
}

// FromTarget implements [Operation]: apply −dt · V.
func (v VelocityGrid) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return v.apply(source, target, lon0, lat0, h0, -v.Dt)
}

func (v VelocityGrid) apply(source, target Spheroid, lon, lat, h, dt float64) (float64, float64, float64, error) {
	if dt == 0 {
		return lon, lat, h, nil
	}
	grid, err := loadVelocityGrid(v.Grid)
	if err != nil {
		return 0, 0, 0, err
	}
	e, n, u, err := grid.Velocity(lon, lat) // mm/year
	if err != nil {
		return 0, 0, 0, err
	}
	// mm/year → m/year
	e *= 0.001
	n *= 0.001
	u *= 0.001

	phi := radian(lat)
	lam := radian(lon)
	sinPhi, cosPhi := math.Sincos(phi)
	sinLam, cosLam := math.Sincos(lam)

	vx := (-sinPhi*cosLam*n - sinLam*e + cosPhi*cosLam*u) * dt
	vy := (-sinPhi*sinLam*n + cosLam*e + cosPhi*sinLam*u) * dt
	vz := (cosPhi*n + sinPhi*u) * dt

	x, y, z := source.GeographicToGeocentric(lon, lat, h)
	lon1, lat1, h1 := target.GeocentricToGeographic(x+vx, y+vy, z+vz)
	return lon1, lat1, h1, nil
}

type velocityGridData struct {
	File    string
	LonUL   float64
	LatUL   float64
	DLon    float64
	DLat    float64
	Columns int
	Rows    int
	East    []float32
	North   []float32
	Up      []float32
}

var velocityGridStore sync.Map

func loadVelocityGrid(name string) (grid *velocityGridData, err error) {
	name = gridFilename(name)

	if v, ok := velocityGridStore.Load(name); ok {
		return v.(*velocityGridData), nil
	}

	defer func() {
		if err == nil {
			grid.File = name
			velocityGridStore.Store(name, grid)
		}
	}()

	r, err := openGridReader(name)
	if err != nil {
		return nil, UnsupportedError{
			Err: fmt.Errorf("velocity grid: %s", name),
		}
	}

	defer r.Close() //nolint:errcheck

	switch strings.ToLower(path.Ext(name)) {
	case ".tif":
		return parseVelocityGridTIFF(r)
	default:
		return nil, fmt.Errorf("crs: unsupported velocity grid extension %q", path.Ext(name))
	}
}

// Velocity returns east/north/up in mm/year at (lon, lat) degrees.
func (g *velocityGridData) Velocity(lon, lat float64) (east, north, up float64, err error) {
	if g.Columns < 2 || g.Rows < 2 || len(g.East) == 0 {
		return 0, 0, 0, UnsupportedError{
			Err: errors.New("invalid grid"),
		}
	}

	colF := (lon - g.LonUL) / g.DLon
	rowF := (g.LatUL - lat) / g.DLat

	col, dx, ok := interpolateIndex(colF, g.Columns)
	if !ok {
		return 0, 0, 0, UnsupportedError{
			Err: errors.New("cannot interpolate grid"),
		}
	}

	row, dy, ok := interpolateIndex(rowF, g.Rows)
	if !ok {
		return 0, 0, 0, UnsupportedError{
			Err: errors.New("cannot interpolate grid"),
		}
	}

	i00 := row*g.Columns + col
	i10 := i00 + 1
	i01 := i00 + g.Columns
	i11 := i01 + 1

	bilinear := func(band []float32) (float64, error) {
		v00 := float64(band[i00])
		v10 := float64(band[i10])
		v01 := float64(band[i01])
		v11 := float64(band[i11])

		if math.IsNaN(v00) || math.IsNaN(v10) || math.IsNaN(v01) || math.IsNaN(v11) {
			return 0, UnsupportedError{
				Err: errors.New("cannot interpolate grid"),
			}
		}

		return (1-dx)*(1-dy)*v00 + dx*(1-dy)*v10 + (1-dx)*dy*v01 + dx*dy*v11, nil
	}

	east, err = bilinear(g.East)
	if err != nil {
		return 0, 0, 0, err
	}

	north, err = bilinear(g.North)
	if err != nil {
		return 0, 0, 0, err
	}

	up, err = bilinear(g.Up)
	if err != nil {
		return 0, 0, 0, err
	}

	return east, north, up, nil
}

func parseVelocityGridTIFF(r io.Reader) (*velocityGridData, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}

	p := &tifParser{data: data}
	if err := p.readHeader(); err != nil {
		return nil, err
	}

	var grids []*velocityGridData
	for ifdOff := p.firstIFD; ifdOff != 0; {
		g, next, err := p.readVelocityIFD(ifdOff)
		if err != nil {
			return nil, err
		}
		if g != nil {
			grids = append(grids, g)
		}
		ifdOff = next
	}

	if len(grids) == 0 {
		return nil, UnsupportedError{
			Err: errors.New("no velocity grid directories"),
		}
	}

	return grids[0], nil
}

func (p *tifParser) readVelocityIFD(off uint32) (*velocityGridData, uint32, error) {
	if int(off)+2 > len(p.data) {
		return nil, 0, fmt.Errorf("tiff: IFD out of range")
	}
	n := int(p.u16(int(off)))
	base := int(off) + 2
	need := base + n*12 + 4
	if need > len(p.data) {
		return nil, 0, fmt.Errorf("tiff: IFD truncated")
	}

	tags := map[uint16]ifdEntry{}
	for i := range n {
		eo := base + i*12
		e := ifdEntry{
			tag:   p.u16(eo),
			typ:   p.u16(eo + 2),
			count: p.u32(eo + 4),
			val:   p.u32(eo + 8),
		}
		tags[e.tag] = e
	}
	next := p.u32(base + n*12)

	meta := parseGDALMetadata(p.tagString(tags, 42112))
	typ := meta["TYPE"]

	switch typ {
	case "VELOCITY", "VELOCITY_CARTOGRAPHIC":
		// ok
	case "", "HORIZONTAL_OFFSET", "VERTICAL_OFFSET_GEOGRAPHIC_TO_VERTICAL", "VERTICAL_OFFSET_VERTICAL_TO_VERTICAL":
		return nil, next, nil // skip non-velocity IFDs
	default:
		return nil, next, fmt.Errorf("tiff: unsupported velocity TYPE %q", typ)
	}

	width := int(p.tagLong(tags, 256))
	height := int(p.tagLong(tags, 257))
	if width <= 0 || height <= 0 {
		return nil, next, fmt.Errorf("tiff: invalid dimensions %dx%d", width, height)
	}

	samples := int(p.tagLong(tags, 277))
	if samples < 3 {
		return nil, next, fmt.Errorf("tiff: velocity grid needs 3 samples, got %d", samples)
	}

	bands, err := p.readBands(tags, width, height, samples)
	if err != nil {
		return nil, next, err
	}

	tie := p.tagDoubles(tags, 33922, 6)
	scale := p.tagDoubles(tags, 33550, 3)
	if len(tie) < 6 || len(scale) < 2 {
		return nil, next, fmt.Errorf("tiff: missing ModelTiepointTag/ModelPixelScaleTag")
	}

	dLon := scale[0]
	dLat := math.Abs(scale[1])
	if dLon == 0 || dLat == 0 {
		return nil, next, fmt.Errorf("tiff: zero pixel scale")
	}

	return &velocityGridData{
		LonUL:   tie[3],
		LatUL:   tie[4],
		DLon:    dLon,
		DLat:    dLat,
		Columns: width,
		Rows:    height,
		East:    bands[0],
		North:   bands[1],
		Up:      bands[2],
	}, next, nil
}
