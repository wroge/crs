package crs

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"path"
	"strings"
	"sync"
)

// VerticalGrid is a geoid / vertical undulation grid (PROJ vgridshift / GTG).
// Grid samples are the geoid undulation N (metres), matching PROJ egm08-style
// grids used with +inv +proj=vgridshift: H = h − N, h = H + N.
type VerticalGrid string

func (vg VerticalGrid) String() string {
	return fmt.Sprintf("operation=vertical_grid grid=%s", string(vg))
}

// ToTarget implements [Operation]: orthometric → ellipsoidal (h = H + N).
func (vg VerticalGrid) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	grid, err := loadVerticalGrid(string(vg))
	if err != nil {
		return 0, 0, 0, err
	}
	n, err := grid.Undulation(lon, lat)
	if err != nil {
		return 0, 0, 0, err
	}
	return lon, lat, h + n, nil
}

// FromTarget implements [Operation]: ellipsoidal → orthometric (H = h − N).
func (vg VerticalGrid) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	grid, err := loadVerticalGrid(string(vg))
	if err != nil {
		return 0, 0, 0, err
	}
	n, err := grid.Undulation(lon0, lat0)
	if err != nil {
		return 0, 0, 0, err
	}
	return lon0, lat0, h0 - n, nil
}

type verticalGridData struct {
	File    string
	LonUL   float64 // degrees east at pixel (0,0) center
	LatUL   float64 // degrees north at pixel (0,0) center
	DLon    float64 // degrees per column (east positive)
	DLat    float64 // degrees per row magnitude (southward in TIFF)
	Columns int
	Rows    int
	Values  []float32 // row-major, TIFF order (row 0 = north)
}

var verticalGridStore sync.Map

func loadVerticalGrid(name string) (grid *verticalGridData, err error) {
	name = gridFilename(name)

	if v, ok := verticalGridStore.Load(name); ok {
		return v.(*verticalGridData), nil
	}

	defer func() {
		if err == nil {
			grid.File = name
			verticalGridStore.Store(name, grid)
		}
	}()

	r, err := openGridReader(name)
	if err != nil {
		return nil, UnsupportedError{
			Err: fmt.Errorf("grid: %s", name),
		}
	}
	defer r.Close() //nolint:errcheck

	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	return parseVerticalGrid(name, data)
}

// selectableOperation reports whether a datum edge can be chosen.
// EPSG still names legacy .grd vertical grids. Those files are not on the
// grid CDN and are not GeoTIFF or GTX, so the next operation is used.
func selectableOperation(op Operation) bool {
	vg, ok := op.(VerticalGrid)
	if !ok {
		return true
	}
	return strings.ToLower(path.Ext(gridFilename(string(vg)))) != ".grd"
}

func parseVerticalGrid(name string, data []byte) (*verticalGridData, error) {
	if looksLikeTIFF(data) {
		return parseVerticalGridTIFF(bytes.NewReader(data))
	}

	switch strings.ToLower(path.Ext(name)) {
	case ".gtx", ".grd":
		return parseVerticalGridGTX(data)
	case ".tif", ".tiff":
		return nil, fmt.Errorf("crs: %s is not a TIFF vertical grid", name)
	default:
		return nil, fmt.Errorf("crs: unsupported vertical grid extension %q", path.Ext(name))
	}
}

func looksLikeTIFF(data []byte) bool {
	if len(data) < 8 {
		return false
	}
	var order binary.ByteOrder
	switch string(data[:2]) {
	case "II":
		order = binary.LittleEndian
	case "MM":
		order = binary.BigEndian
	default:
		return false
	}
	return order.Uint16(data[2:4]) == 42
}

// parseVerticalGridGTX reads a NOAA GTX geoid grid (big-endian).
// Row 0 of the file is the southern edge. The returned grid uses the same
// north-first layout as a GeoTIFF vertical grid.
func parseVerticalGridGTX(data []byte) (*verticalGridData, error) {
	if len(data) < 40 {
		return nil, fmt.Errorf("gtx: short header")
	}

	yOrigin := math.Float64frombits(binary.BigEndian.Uint64(data[0:8]))
	xOrigin := math.Float64frombits(binary.BigEndian.Uint64(data[8:16]))
	yStep := math.Float64frombits(binary.BigEndian.Uint64(data[16:24]))
	xStep := math.Float64frombits(binary.BigEndian.Uint64(data[24:32]))
	rows := int(int32(binary.BigEndian.Uint32(data[32:36])))
	cols := int(int32(binary.BigEndian.Uint32(data[36:40])))

	if cols <= 0 || rows <= 0 || xStep <= 0 || yStep <= 0 ||
		xOrigin < -360 || xOrigin > 360 || yOrigin < -90 || yOrigin > 90 {
		return nil, fmt.Errorf("gtx: invalid header")
	}
	if xOrigin >= 180 {
		xOrigin -= 360
	}

	if cols > (len(data)-40)/4/rows {
		return nil, fmt.Errorf("gtx: truncated grid")
	}
	n := rows * cols

	values := make([]float32, n)
	for row := range rows {
		src := (rows - 1 - row) * cols
		dst := row * cols
		for col := range cols {
			off := 40 + (src+col)*4
			v := math.Float32frombits(binary.BigEndian.Uint32(data[off : off+4]))
			if gtxNodata(v) {
				v = float32(math.NaN())
			}
			values[dst+col] = v
		}
	}

	return &verticalGridData{
		LonUL:   xOrigin,
		LatUL:   yOrigin + float64(rows-1)*yStep,
		DLon:    xStep,
		DLat:    yStep,
		Columns: cols,
		Rows:    rows,
		Values:  values,
	}, nil
}

func gtxNodata(v float32) bool {
	return v > 1000 || v < -1000 || v == float32(-88.88880)
}

// Undulation returns geoid undulation N in metres at (lon, lat) degrees.
func (g *verticalGridData) Undulation(lon, lat float64) (float64, error) {
	if g.Columns < 2 || g.Rows < 2 || len(g.Values) == 0 {
		return 0, UnsupportedError{
			Err: errors.New("invalid grid"),
		}
	}

	colF := (lon - g.LonUL) / g.DLon
	rowF := (g.LatUL - lat) / g.DLat

	col, dx, ok := interpolateIndex(colF, g.Columns)
	if !ok {
		return 0, UnsupportedError{
			Err: errors.New("cannot interpolate grid"),
		}
	}
	row, dy, ok := interpolateIndex(rowF, g.Rows)
	if !ok {
		return 0, UnsupportedError{
			Err: errors.New("cannot interpolate grid"),
		}
	}

	i00 := row*g.Columns + col
	i10 := i00 + 1
	i01 := i00 + g.Columns
	i11 := i01 + 1

	v00 := float64(g.Values[i00])
	v10 := float64(g.Values[i10])
	v01 := float64(g.Values[i01])
	v11 := float64(g.Values[i11])

	if math.IsNaN(v00) || math.IsNaN(v10) || math.IsNaN(v01) || math.IsNaN(v11) {
		return 0, UnsupportedError{
			Err: errors.New("cannot interpolate grid"),
		}
	}

	n := (1-dx)*(1-dy)*v00 + dx*(1-dy)*v10 + (1-dx)*dy*v01 + dx*dy*v11
	return n, nil
}

// parseVerticalGridTIFF reads a PROJ vertical/geoid GeoTIFF (GTG).
func parseVerticalGridTIFF(r io.Reader) (*verticalGridData, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}

	p := &tifParser{data: data}
	if err := p.readHeader(); err != nil {
		return nil, err
	}

	var grids []*verticalGridData

	for ifdOff := p.firstIFD; ifdOff != 0; {
		g, next, err := p.readVerticalIFD(ifdOff)
		if err != nil {
			return nil, err
		}
		if g != nil {
			grids = append(grids, g)
		}
		ifdOff = next
	}

	if len(grids) == 0 {
		return nil, fmt.Errorf("tiff: no vertical grid directories")
	}

	return grids[0], nil
}

func (p *tifParser) readVerticalIFD(off uint32) (*verticalGridData, uint32, error) {
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
	case "", "VERTICAL_OFFSET_GEOGRAPHIC_TO_VERTICAL", "VERTICAL_OFFSET_VERTICAL_TO_VERTICAL":
		// ok
	case "HORIZONTAL_OFFSET":
		return nil, next, nil // skip horizontal IFDs in mixed files
	default:
		return nil, next, fmt.Errorf("tiff: unsupported vertical TYPE %q", typ)
	}

	width := int(p.tagLong(tags, 256))
	height := int(p.tagLong(tags, 257))
	if width <= 0 || height <= 0 {
		return nil, next, fmt.Errorf("tiff: invalid dimensions %dx%d", width, height)
	}

	samples := max(int(p.tagLong(tags, 277)), 1)

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

	return &verticalGridData{
		LonUL:   tie[3],
		LatUL:   tie[4],
		DLon:    dLon,
		DLat:    dLat,
		Columns: width,
		Rows:    height,
		Values:  bands[0],
	}, next, nil
}
