package crs

import (
	"encoding/binary"
	"math"
	"testing"
)

func TestParseVerticalGridGTX(t *testing.T) {
	// 2×2 nodes. File row 0 is south: 1, 2 at lat 10; 3, 4 at lat 11.
	data := gtxGrid(10, 20, 1, 1, []float32{1, 2, 3, 4})
	g, err := parseVerticalGridGTX(data)
	if err != nil {
		t.Fatal(err)
	}
	if g.LatUL != 11 || g.LonUL != 20 || g.Rows != 2 || g.Columns != 2 {
		t.Fatalf("header lat %v lon %v %dx%d", g.LatUL, g.LonUL, g.Columns, g.Rows)
	}
	n, err := g.Undulation(20.5, 10.5)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(n-2.5) > 1e-6 {
		t.Fatalf("undulation %v", n)
	}

	data[40] = 0
	binary.BigEndian.PutUint32(data[40:44], math.Float32bits(float32(-88.88880)))
	g, err = parseVerticalGridGTX(data)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := g.Undulation(20, 10); err == nil {
		t.Fatal("nodata should not interpolate")
	}
}

func TestSelectableOperationSkipsGRD(t *testing.T) {
	if selectableOperation(VerticalGrid("geo_igm_mar06.grd")) {
		t.Fatal(".grd must not be selected")
	}
	if !selectableOperation(VerticalGrid("us_nga_egm08_25.tif")) {
		t.Fatal(".tif must stay selectable")
	}
	if !selectableOperation(VerticalGrid("Und_min1x1_egm2008_isw=82_WGS84_TideFree")) {
		t.Fatal("extensionless grids must stay selectable")
	}
	if !selectableOperation(VerticalOffset{Dh: 0.141}) {
		t.Fatal("vertical offset must stay selectable")
	}
}

func gtxGrid(yOrigin, xOrigin, yStep, xStep float64, southFirst []float32) []byte {
	cols := 2
	rows := len(southFirst) / cols
	buf := make([]byte, 40+len(southFirst)*4)
	binary.BigEndian.PutUint64(buf[0:8], math.Float64bits(yOrigin))
	binary.BigEndian.PutUint64(buf[8:16], math.Float64bits(xOrigin))
	binary.BigEndian.PutUint64(buf[16:24], math.Float64bits(yStep))
	binary.BigEndian.PutUint64(buf[24:32], math.Float64bits(xStep))
	binary.BigEndian.PutUint32(buf[32:36], uint32(rows))
	binary.BigEndian.PutUint32(buf[36:40], uint32(cols))
	for i, v := range southFirst {
		binary.BigEndian.PutUint32(buf[40+i*4:], math.Float32bits(v))
	}
	return buf
}
