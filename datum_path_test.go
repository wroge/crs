package crs

import (
	"fmt"
	"strings"
	"testing"
)

// TestDatumPaths locks the operation chosen for the grid and vertical cases
// in the PROJ goldens. A later change of accuracy or bbox must show up here
// before it quietly swaps a grid.
func TestDatumPaths(t *testing.T) {
	cases := []struct {
		name            string
		from, to        int
		lon, lat, epoch float64
		want            string
	}{
		{
			name: "horizontal_grid",
			from: 4326, to: 27200,
			lon: 173.70813241, lat: -42.06287211,
			want: "wgs84 -> nzgd49 (inverse operation=horizontal_grid grid=nz_linz_nzgd2kgrid0005.tif)",
		},
		{
			name: "vertical_grid",
			from: 4979, to: 9518,
			lon: -21.0815164, lat: -22.06923813,
			want: "wgs84 -> wgs84_egm2008geoid (inverse operation=vertical_grid grid=Und_min1x1_egm2008_isw=82_WGS84_TideFree)",
		},
		{
			name: "vertical_offset",
			from: 4937, to: 9724,
			lon: 13.63490925, lat: 37.5206217,
			want: "etrs89 -> etrs89_catania1965 (inverse operation=vertical_offset dh=0.141)",
		},
		{
			name: "velocity_grid",
			from: 9067, to: 4669,
			lon: 21.83319205, lat: 55.66216944, epoch: 2010,
			want: "etrf2000 -> etrs89 (identity) | etrs89 -> etrf2014 (inverse identity) | etrf2014 -> nkg_etrf14 (velocity_grid grid=eur_nkg_nkgrf17vel.tif epoch=2000) | nkg_etrf14 -> lks94 (position_vector tx=0.36749 ty=0.14351 tz=-0.18472)",
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, err := datumHops(c.from, c.to, c.lon, c.lat, c.epoch)
			if err != nil {
				t.Fatal(err)
			}
			if c.want == "" {
				t.Fatalf("lock this path:\n%s", got)
			}
			if got != c.want {
				t.Fatalf("path changed\n got %s\nwant %s", got, c.want)
			}
		})
	}
}

func datumHops(fromCode, toCode int, lon, lat, epoch float64) (string, error) {
	from, err := Load(fromCode)
	if err != nil {
		return "", err
	}
	to, err := Load(toCode)
	if err != nil {
		return "", err
	}
	hasEpoch := epoch != 0
	path, err := findBestPath(from.Datum, to.Datum, lon, lat, nil, hasEpoch, epoch)
	if err != nil {
		return "", err
	}
	if len(path) < 2 {
		return strings.Join(path, " = "), nil
	}

	var parts []string
	for i := 0; i < len(path)-1; i++ {
		cands := edgesBetween(path[i], path[i+1], lon, lat, nil, hasEpoch, epoch)
		if len(cands) == 0 {
			parts = append(parts, path[i]+" = "+path[i+1])
			continue
		}
		e := cands[0]
		op := "?"
		owner, err := loadDatum(e.key.owner)
		if err == nil && e.key.index >= 0 && e.key.index < len(owner.Transformations) {
			op = fmt.Sprint(owner.Transformations[e.key.index].Operation)
			if e.inverse {
				op = "inverse " + op
			}
		}
		parts = append(parts, path[i]+" -> "+path[i+1]+" ("+op+")")
	}
	return strings.Join(parts, " | "), nil
}
