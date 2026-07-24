package crs

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"strings"
)

func parseGridShiftBinary(r io.Reader) (*horizontalGridData, error) {
	first, err := readGridRecord(r)
	if err != nil {
		return nil, fmt.Errorf("read file header: %w", err)
	}

	order := detectGSBEndian(first)

	numSrec, numFile, err := parseGridFileHeader(first, r, order)
	if err != nil {
		return nil, fmt.Errorf("read file header: %w", err)
	}

	n := []subGrid{}

	for range numFile {
		sg, count, err := readGridSubgridHeader(r, numSrec, order)
		if err != nil {
			return nil, fmt.Errorf("read subgrid header: %w", err)
		}

		sg.Values = make([][2]float32, count)

		for i := range sg.Values {
			rec, err := readGridRecord(r)
			if err != nil {
				return nil, fmt.Errorf("read grid value %d: %w", i, err)
			}

			sg.Values[i][0] = recordFloat32(rec, 0, order)
			sg.Values[i][1] = recordFloat32(rec, 4, order)
		}

		if err := sg.validate(); err != nil {
			return nil, err
		}

		n = append(n, sg)
	}

	return &horizontalGridData{
		subGrids: n,
	}, nil
}

type gridRecord [16]byte

func detectGSBEndian(rec gridRecord) binary.ByteOrder {
	if recordName(rec) != "NUM_OREC" {
		return binary.LittleEndian
	}

	le := int32(binary.LittleEndian.Uint32(rec[8:12]))
	if le > 0 && le < 100 {
		return binary.LittleEndian
	}

	be := int32(binary.BigEndian.Uint32(rec[8:12]))
	if be > 0 && be < 100 {
		return binary.BigEndian
	}

	return binary.LittleEndian
}

func parseGridFileHeader(first gridRecord, r io.Reader, order binary.ByteOrder) (numSrec, numFile int, err error) {
	rec := first
	numOrec := int(recordValueInt32(rec, order))

	if recordName(rec) == "NUM_OREC" && numOrec > 0 {
		for i := 1; i < numOrec; i++ {
			rec, err = readGridRecord(r)
			if err != nil {
				return 0, 0, err
			}

			switch recordName(rec) {
			case "NUM_SREC":
				numSrec = int(recordValueInt32(rec, order))
			case "NUM_FILE":
				numFile = int(recordValueInt32(rec, order))
			}
		}
	}

	if numSrec <= 0 {
		numSrec = 11
	}

	return numSrec, numFile, nil
}

func readGridSubgridHeader(r io.Reader, numSrec int, order binary.ByteOrder) (subGrid, int, error) {
	var (
		sg    subGrid
		count int
		cn    [4]string
	)

	for range numSrec {
		rec, err := readGridRecord(r)
		if err != nil {
			return subGrid{}, 0, err
		}

		switch recordName(rec) {
		case "SUB_NAME":
			sg.Name = recordValueString(rec)
		case "CN_1":
			cn[0] = recordValueString(rec)
		case "CN_2":
			cn[1] = recordValueString(rec)
		case "CN_3":
			cn[2] = recordValueString(rec)
		case "CN_4":
			cn[3] = recordValueString(rec)
		case "PARENT":
			sg.Parent = recordValueString(rec)
		case "S_LAT":
			sg.SLat = recordValueFloat64(rec, order)
		case "N_LAT":
			sg.NLat = recordValueFloat64(rec, order)
		case "E_LONG":
			sg.ELong = recordValueFloat64(rec, order)
		case "W_LONG":
			sg.WLong = recordValueFloat64(rec, order)
		case "LAT_INC":
			sg.LatInc = recordValueFloat64(rec, order)
		case "LONG_INC":
			sg.LongInc = recordValueFloat64(rec, order)
		case "COLUMNS":
			sg.Columns = int(recordValueInt32(rec, order))
		case "ROWS":
			sg.Rows = int(recordValueInt32(rec, order))
		case "GS_COUNT":
			count = int(recordValueInt32(rec, order))
		}
	}

	if sg.Name == "" {
		sg.Name = trimGrid(strings.Join(cn[:], ""))
	}

	if count <= 0 {
		return subGrid{}, 0, fmt.Errorf("subgrid %q: missing GS_COUNT", sg.Name)
	}

	if sg.Columns == 0 || sg.Rows == 0 {
		if sg.LatInc == 0 || sg.LongInc == 0 {
			return subGrid{}, 0, fmt.Errorf("subgrid %q: missing grid dimensions", sg.Name)
		}

		sg.Rows = int(math.Floor((sg.NLat-sg.SLat)/sg.LatInc+0.5)) + 1
		sg.Columns = int(math.Floor((sg.WLong-sg.ELong)/sg.LongInc+0.5)) + 1
	}

	if sg.Columns*sg.Rows != count {
		return subGrid{}, 0, fmt.Errorf("subgrid %q: GS_COUNT %d != %d×%d", sg.Name, count, sg.Columns, sg.Rows)
	}

	return sg, count, nil
}

func trimGrid(s string) string {
	return strings.TrimFunc(s, func(r rune) bool {
		return r == 0 || r == ' '
	})
}

func recordValueString(rec gridRecord) string {
	return trimGrid(string(rec[8:16]))
}

func recordValueInt32(rec gridRecord, order binary.ByteOrder) int32 {
	return int32(order.Uint32(rec[8:12]))
}

func recordValueFloat64(rec gridRecord, order binary.ByteOrder) float64 {
	return math.Float64frombits(order.Uint64(rec[8:16]))
}

func readGridRecord(r io.Reader) (gridRecord, error) {
	var rec gridRecord

	_, err := io.ReadFull(r, rec[:])

	return rec, err
}

func recordName(rec gridRecord) string {
	return trimGrid(string(rec[:8]))
}

func recordFloat32(rec gridRecord, off int, order binary.ByteOrder) float32 {
	return math.Float32frombits(order.Uint32(rec[off : off+4]))
}
