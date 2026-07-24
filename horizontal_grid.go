package crs

import (
	"fmt"
)

type HorizontalGrid string

func (hg HorizontalGrid) String() string {
	return fmt.Sprintf("operation=horizontal_grid grid=%s", string(hg))
}

func (hg HorizontalGrid) FromTarget(source Spheroid, target Spheroid, lon0 float64, lat0 float64, h0 float64) (float64, float64, float64, error) {
	grid, err := loadGrid(string(hg))
	if err != nil {
		return 0, 0, 0, err
	}

	return grid.FromTarget(source, target, lon0, lat0, h0)
}

func (hg HorizontalGrid) ToTarget(source Spheroid, target Spheroid, lon float64, lat float64, h float64) (float64, float64, float64, error) {
	grid, err := loadGrid(string(hg))
	if err != nil {
		return 0, 0, 0, err
	}

	return grid.ToTarget(source, target, lon, lat, h)
}

type horizontalGridData struct {
	File     string
	subGrids []subGrid
}

func (g *horizontalGridData) FromTarget(source Spheroid, target Spheroid, lon0 float64, lat0 float64, h0 float64) (float64, float64, float64, error) {
	lon, lat, err := g.FromWGS84(lon0, lat0)
	if err != nil {
		return 0, 0, 0, err
	}

	return lon, lat, h0, nil
}

func (g *horizontalGridData) ToTarget(source Spheroid, target Spheroid, lon float64, lat float64, h float64) (float64, float64, float64, error) {
	lon0, lat0, err := g.ToWGS84(lon, lat)
	if err != nil {
		return 0, 0, 0, err
	}

	return lon0, lat0, h, nil
}

type Inverse struct {
	Operation Operation
}

func (i Inverse) String() string {
	return fmt.Sprintf("%s inverse", i.Operation)
}

func (i Inverse) FromTarget(source Spheroid, target Spheroid, lon0 float64, lat0 float64, h0 float64) (float64, float64, float64, error) {
	return i.Operation.ToTarget(target, source, lon0, lat0, h0)
}

func (i Inverse) ToTarget(source Spheroid, target Spheroid, lon float64, lat float64, h float64) (float64, float64, float64, error) {
	return i.Operation.FromTarget(target, source, lon, lat, h)
}
