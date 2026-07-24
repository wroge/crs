# crs

[![Coverage](https://wroge.github.io/crs/badges/coverage.svg)](https://wroge.github.io/crs/)
[![go.dev reference](https://img.shields.io/badge/go.dev-reference-007d9c?logo=go&logoColor=white)](https://pkg.go.dev/github.com/wroge/crs@v0.0.2)

```sh
go get github.com/wroge/crs@v0.0.2
```

Go library for coordinate reference systems, map projections, and datum
transformations. Supports the majority of EPSG codes (thousands of embedded
definitions).

```go
transform, err := crs.Transform(4326, 32632)
if err != nil {
	log.Fatal(err)
}

east, north, h, err := transform(10.0, 50.0, 0)
```

Axis conventions: geographic `lon, lat, h`; projected `east, north, h`;
geocentric `x, y, z`.

## Geodetic grids

Some datum (nad83, osgb36, ...) transforms need grid files. Provide them via CDN and/or local sources:

```go
crs.SetGridCDN("https://cdn.proj.org/")                          // download on demand
crs.SetGridCDN("https://cdn.proj.org/", "./.cache") // optional disk cache

crs.RegisterGridDir("/path/to/grids") // search a directory
crs.RegisterGridFS("", embedFS)       // or an fs.FS (e.g. embed.FS)
```

Lookup order: registered directories, then filesystems, then CDN.

## Time-dependent transforms

```go
transform, err := crs.TransformAt(9989, 9069, 2010.0) // coordinate epoch (decimal year)
```

## License

Licensed under the [MIT License](LICENSE).
