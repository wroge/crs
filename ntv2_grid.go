package crs

import (
	"bytes"
	"fmt"
	"io"
	"io/fs"
	"math"
	"net/http"
	"os"
	"path"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"sync"
	"time"
)

type gridFilesystem struct {
	prefix string
	fsys   fs.FS
}

var (
	gridMu         sync.RWMutex
	gridCDNBase    string
	gridCDNCache   string
	gridSearchDirs []string
	gridFS         []gridFilesystem // sorted by prefix length, longest first
	gridHTTPClient = &http.Client{Timeout: 10 * time.Minute}
)

func SetGridCDN(uri string, cacheDir ...string) {
	gridMu.Lock()
	defer gridMu.Unlock()

	uri = strings.TrimSpace(uri)
	if uri == "" {
		gridCDNBase = ""
		gridCDNCache = ""
		return
	}

	if !strings.HasSuffix(uri, "/") {
		uri += "/"
	}

	gridCDNBase = uri

	gridCDNCache = ""
	if len(cacheDir) > 0 {
		gridCDNCache = strings.TrimSpace(cacheDir[0])
	}
}

func RegisterGridDir(dir string) {
	dir = strings.TrimSpace(dir)
	if dir == "" {
		return
	}

	gridMu.Lock()
	gridSearchDirs = append(gridSearchDirs, dir)
	gridMu.Unlock()
}

func RegisterGridFS(prefix string, fsys fs.FS) {
	if fsys == nil {
		panic("crs: RegisterGridFS: fsys is nil")
	}

	if prefix == "/" {
		prefix = ""
	}

	gridMu.Lock()

	gridFS = append(gridFS, gridFilesystem{prefix: prefix, fsys: fsys})

	sort.Slice(gridFS, func(i, j int) bool {
		return len(gridFS[i].prefix) > len(gridFS[j].prefix)
	})

	gridMu.Unlock()
}

func gridFilename(p string) string {
	return path.Base(strings.ReplaceAll(p, `\`, `/`))
}

func openGridReader(name string) (io.ReadCloser, error) {
	gridMu.RLock()

	cdnBase := gridCDNBase
	cacheDir := gridCDNCache
	searchDirs := append([]string(nil), gridSearchDirs...)
	filesystems := append([]gridFilesystem(nil), gridFS...)

	gridMu.RUnlock()

	if cacheDir != "" {
		found := slices.Contains(searchDirs, cacheDir)
		if !found {
			searchDirs = append(searchDirs, cacheDir)
		}
	}

	for _, dir := range searchDirs {
		if dir == "" {
			continue
		}

		if f, err := os.Open(filepath.Join(dir, name)); err == nil {
			return f, nil
		}
	}

	for _, reg := range filesystems {
		p := name
		if reg.prefix != "" {
			p = strings.TrimSuffix(reg.prefix, "/") + "/" + name
		}

		if f, err := reg.fsys.Open(p); err == nil {
			return f, nil
		}
	}

	if cdnBase == "" {
		return nil, fmt.Errorf("crs: grid %q not found", name)
	}

	resp, err := gridHTTPClient.Get(cdnBase + name)
	if err != nil {
		return nil, fmt.Errorf("crs: grid %q not found", name)
	}

	if resp.StatusCode == http.StatusNotFound {
		_ = resp.Body.Close()
		return nil, fmt.Errorf("crs: grid %q not found", name)
	}

	if resp.StatusCode != http.StatusOK {
		_ = resp.Body.Close()
		return nil, fmt.Errorf("cdn: grid %q: HTTP %s", name, resp.Status)
	}

	if cacheDir == "" {
		return resp.Body, nil
	}

	data, err := io.ReadAll(resp.Body)
	_ = resp.Body.Close()
	if err != nil {
		return nil, err
	}

	final := filepath.Join(cacheDir, name)
	tmp := final + ".tmp"

	if err = os.MkdirAll(cacheDir, 0o755); err != nil {
		return io.NopCloser(bytes.NewReader(data)), nil
	}

	if err = os.WriteFile(tmp, data, 0o644); err != nil {
		return io.NopCloser(bytes.NewReader(data)), nil
	}

	if err = os.Rename(tmp, final); err != nil {
		_ = os.Remove(tmp)

		return io.NopCloser(bytes.NewReader(data)), nil
	}

	if f, err := os.Open(final); err == nil {
		return f, nil
	}

	return io.NopCloser(bytes.NewReader(data)), nil
}

// gridStore caches parsed grids in memory (local embed, disk, and CDN sources).
var gridStore sync.Map

func loadGrid(name string) (grid *horizontalGridData, err error) {
	name = gridFilename(name)

	if v, ok := gridStore.Load(name); ok {
		return v.(*horizontalGridData), nil
	}

	defer func() {
		if err == nil {
			grid.File = name
			gridStore.Store(name, grid)
		}
	}()

	r, err := openGridReader(name)
	if err != nil {
		return nil, GridNotFoundError{
			Err: fmt.Errorf("grid: %s", name),
		}
	}

	defer r.Close() //nolint:errcheck

	switch strings.ToLower(path.Ext(name)) {
	case ".tif":
		return parseTaggedImageFileFormat(r)
	case ".gsb":
		return parseGridShiftBinary(r)
	default:
		return nil, fmt.Errorf("crs: unknown grid extension %q", path.Ext(name))
	}
}

type subGrid struct {
	Name    string
	Parent  string
	Columns int
	Rows    int
	SLat    float64
	NLat    float64
	ELong   float64
	WLong   float64
	LatInc  float64
	LongInc float64
	Values  [][2]float32
}

func (g *horizontalGridData) ToWGS84(lon, lat float64) (float64, float64, error) {
	dlon, dlat, err := g.Shift(lon, lat)
	if err != nil {
		return 0, 0, err
	}

	return lon + dlon, lat + dlat, nil
}

func (g *horizontalGridData) FromWGS84(lon, lat float64) (float64, float64, error) {
	qlon, qlat := lon, lat

	for range 10 {
		dlon, dlat, err := g.Shift(qlon, qlat)
		if err != nil {
			return 0, 0, err
		}

		newLon := lon - dlon
		newLat := lat - dlat

		if math.Abs(newLon-qlon) < 1e-12 && math.Abs(newLat-qlat) < 1e-12 {
			break
		}

		qlon, qlat = newLon, newLat
	}

	return qlon, qlat, nil
}

func (g *horizontalGridData) Shift(lon, lat float64) (dlon, dlat float64, err error) {
	idx := g.selectSubgrid(lon, lat)
	if idx < 0 {
		return 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("coordinate: [%f, %f]: %s", lon, lat, g.File),
		}
	}

	return g.subGrids[idx].shift(lon, lat)
}

func (g *horizontalGridData) selectSubgrid(lon, lat float64) int {
	lam := -lon * 3600
	phi := lat * 3600

	for i, sg := range g.subGrids {
		if sg.Parent != "" && sg.Parent != "NONE" {
			continue
		}
		if sg.contains(phi, lam) {
			return g.deepestSubgridRec(i, phi, lam, make(map[int]bool))
		}
	}

	return -1
}

func (g *horizontalGridData) deepestSubgridRec(idx int, phi, lam float64, seen map[int]bool) int {
	if seen[idx] {
		return idx
	}

	seen[idx] = true

	best := idx
	name := g.subGrids[idx].Name

	for i, sg := range g.subGrids {
		if sg.Parent != name || !sg.contains(phi, lam) {
			continue
		}

		if deeper := g.deepestSubgridRec(i, phi, lam, seen); deeper >= 0 {
			best = deeper
		}
	}

	return best
}

const relToleranceHGridShift = 1e-5

func (sg subGrid) gridEpsilon() float64 {
	return (sg.LongInc + sg.LatInc) * relToleranceHGridShift
}

func (sg subGrid) contains(phi, lam float64) bool {
	eps := sg.gridEpsilon()

	return phi >= sg.SLat-eps && phi <= sg.NLat+eps &&
		lam >= sg.ELong-eps && lam <= sg.WLong+eps
}

func interpolateIndex(f float64, size int) (idx int, frac float64, ok bool) {
	if math.IsNaN(f) {
		return 0, 0, false
	}

	idx = int(math.Round(math.Floor(f)))
	frac = f - float64(idx)

	if idx < 0 {
		if idx == -1 && frac > 1-10*relToleranceHGridShift {
			idx++
			frac = 0
		} else {
			return 0, 0, false
		}
	} else if idx+1 >= size {
		if idx+1 == size && frac < 10*relToleranceHGridShift {
			idx--
			frac = 1
		} else {
			return 0, 0, false
		}
	}

	return idx, frac, true
}

func (sg subGrid) shift(lon, lat float64) (dlon, dlat float64, err error) {
	if sg.Columns < 2 || sg.Rows < 2 || len(sg.Values) == 0 {
		return 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("coordinate: [%f, %f]: %s", lon, lat, sg.Name),
		}
	}

	phi := lat * 3600
	lam := -lon * 3600
	if !sg.contains(phi, lam) {
		return 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("coordinate: [%f, %f]: %s", lon, lat, sg.Name),
		}
	}

	fcol := (lam - sg.ELong) / sg.LongInc
	frow := (phi - sg.SLat) / sg.LatInc

	ppr := sg.Columns

	col, dx, ok := interpolateIndex(fcol, ppr)
	if !ok {
		return 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("coordinate: [%f, %f]: %s", lon, lat, sg.Name),
		}
	}

	row, dy, ok := interpolateIndex(frow, sg.Rows)
	if !ok {
		return 0, 0, OutOfBoundsError{
			Err: fmt.Errorf("coordinate: [%f, %f]: %s", lon, lat, sg.Name),
		}
	}

	se := row*ppr + col
	sw := se + 1
	ne := se + ppr
	nw := ne + 1

	sse := sg.Values[se]
	ssw := sg.Values[sw]
	sne := sg.Values[ne]
	snw := sg.Values[nw]

	latsv := (1-dx)*(1-dy)*float64(sse[0]) + dx*(1-dy)*float64(ssw[0]) +
		(1-dx)*dy*float64(sne[0]) + dx*dy*float64(snw[0])
	lonsv := (1-dx)*(1-dy)*float64(sse[1]) + dx*(1-dy)*float64(ssw[1]) +
		(1-dx)*dy*float64(sne[1]) + dx*dy*float64(snw[1])

	return -lonsv / 3600, latsv / 3600, nil
}

func (sg subGrid) validate() error {
	want := sg.Columns * sg.Rows
	if want == 0 {
		return fmt.Errorf("subgrid %q: zero grid dimensions", sg.Name)
	}

	if len(sg.Values) != want {
		return fmt.Errorf("subgrid %q: got %d values, want %d", sg.Name, len(sg.Values), want)
	}

	return nil
}
