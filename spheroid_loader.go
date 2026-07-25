package crs

import (
	"embed"
	"fmt"
	"io"
	"strings"
	"sync"
)

//go:embed spheroid/*.txt
var spheroidDir embed.FS

var spheroidStore sync.Map

func RegisterSpheroid(s Spheroid) {
	spheroidStore.Store(s.Name, s)
}

func loadSpheroid(name string) (s Spheroid, err error) {
	name = strings.ToLower(strings.TrimSpace(name))

	v, ok := spheroidStore.Load(name)
	if ok {
		return v.(Spheroid), nil
	}

	defer func() {
		if err == nil {
			spheroidStore.Store(name, s)
		}
	}()

	file, err := spheroidDir.Open(fmt.Sprintf("spheroid/%s.txt", name))
	if err != nil {
		return s, err
	}

	data, err := io.ReadAll(file)
	if err != nil {
		return s, err
	}

	s, err = parseSpheroid(string(data))
	if err != nil {
		return s, err
	}

	s.Name = name

	return s, nil
}
