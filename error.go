package crs

import (
	"fmt"
)

type UnsupportedError struct {
	Err error
}

func (e UnsupportedError) Unwrap() error {
	return e.Err
}

func (e UnsupportedError) Error() string {
	return fmt.Sprintf("unsupported: %s", e.Err)
}

type OutOfBoundsError struct {
	Err error
}

func (e OutOfBoundsError) Unwrap() error {
	return e.Err
}

func (e OutOfBoundsError) Error() string {
	return fmt.Sprintf("out of bounds: %s", e.Err)
}

type GridNotFoundError struct {
	Err error
}

func (e GridNotFoundError) Unwrap() error {
	return e.Err
}

func (e GridNotFoundError) Error() string {
	return fmt.Sprintf("grid not found: %s", e.Err)
}
