# crs

[![Coverage](https://wroge.github.io/crs/badges/coverage.svg)](https://wroge.github.io/crs/)

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

## Time-dependent transforms

```go
transform, err := crs.TransformAt(9989, 9069, 2010.0) // coordinate epoch (decimal year)
```

## Geodetic grids

Some datum transforms need grid files. Provide them via CDN and/or local sources:

```go
crs.SetGridCDN("https://cdn.proj.org/")                          // download on demand
crs.SetGridCDN("https://cdn.proj.org/", "./.cache") // optional disk cache

crs.RegisterGridDir("/path/to/grids") // search a directory
crs.RegisterGridFS("", embedFS)       // or an fs.FS (e.g. embed.FS)
```

Lookup order: registered directories, then filesystems, then CDN.

## License

Licensed under the [MIT License](LICENSE).
