package crs

import (
	"container/heap"
	"errors"
	"fmt"
	"io/fs"
	"strings"
	"sync"
)

const maxPathHops = 8

// isHubEquivalent reports whether d may be treated as a free identity hop to/from
// WGS84. Datums with no exported ops (typically PROJ-ballpark-only in the wild)
// stay hub-equivalent so CRS remain reachable without exporting ballpark edges.
// Published null hubs with accuracy/bbox use operation=identity instead of empty.
func isHubEquivalent(d Datum) bool {
	return len(d.Transformations) == 0
}

type revEdge struct {
	owner string
	index int
}

type edgeKey struct {
	owner   string
	index   int
	inverse bool
}

type graphEdge struct {
	key      edgeKey
	to       string
	accuracy float64
	inverse  bool
}

var (
	reverseOnce  sync.Once
	reverseIndex map[string][]revEdge
)

func ensureReverseIndex() {
	reverseOnce.Do(func() {
		reverseIndex = make(map[string][]revEdge)

		entries, err := fs.ReadDir(datumDir, "datum")
		if err != nil {
			return
		}

		for _, e := range entries {
			if e.IsDir() || !strings.HasSuffix(e.Name(), ".txt") {
				continue
			}

			name := strings.TrimSuffix(e.Name(), ".txt")
			d, err := loadDatum(name)
			if err != nil {
				continue
			}

			for i, t := range d.Transformations {
				if t.Target == nil || t.Operation == nil {
					continue
				}

				tgt := strings.ToLower(t.Target.Name)
				reverseIndex[tgt] = append(reverseIndex[tgt], revEdge{owner: name, index: i})
			}
		}
	})
}

func datumName(d Datum) string {
	return strings.ToLower(d.Name)
}

// outgoingEdges returns graph edges from node that cover (lon, lat), excluding keys in excluded.
// Forward and inverse candidates are already ordered by the owner's transformation sort
// (accuracy, then area); we emit one best edge per neighbor for Dijkstra, but
// edgesBetween returns all candidates for apply.
// Time-specific Helmerts are included only when hasCoordEpoch and coordEpoch match.
func outgoingEdges(node string, lon, lat float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) []graphEdge {
	ensureReverseIndex()

	d, err := loadDatum(node)
	if err != nil {
		return nil
	}

	// AoU bboxes are Greenwich-based; lon may be relative to this node's PM.
	glon := lon + d.primeMeridianLongitude()

	best := make(map[string]graphEdge) // neighbor -> best edge

	consider := func(e graphEdge) {
		if excluded[e.key] {
			return
		}
		cur, ok := best[e.to]
		if !ok || e.accuracy < cur.accuracy {
			best[e.to] = e
		}
	}

	for i, t := range d.Transformations {
		if t.Target == nil || t.Operation == nil {
			continue
		}

		if !timeSpecificAllowed(t.Operation, hasCoordEpoch, coordEpoch) {
			continue
		}

		if !t.BoundingBox.Contains(glon, lat) {
			continue
		}

		tgt := strings.ToLower(t.Target.Name)
		if tgt == "" || tgt == node {
			continue
		}

		consider(graphEdge{
			key:      edgeKey{owner: node, index: i, inverse: false},
			to:       tgt,
			accuracy: t.Accuracy,
			inverse:  false,
		})
	}

	for _, rev := range reverseIndex[node] {
		owner, err := loadDatum(rev.owner)
		if err != nil || rev.index < 0 || rev.index >= len(owner.Transformations) {
			continue
		}

		t := owner.Transformations[rev.index]
		if t.Target == nil || t.Operation == nil {
			continue
		}

		if !timeSpecificAllowed(t.Operation, hasCoordEpoch, coordEpoch) {
			continue
		}

		if !t.BoundingBox.Contains(glon, lat) {
			continue
		}

		consider(graphEdge{
			key:      edgeKey{owner: rev.owner, index: rev.index, inverse: true},
			to:       rev.owner,
			accuracy: t.Accuracy,
			inverse:  true,
		})
	}

	out := make([]graphEdge, 0, len(best))
	for _, e := range best {
		out = append(out, e)
	}

	return out
}

// edgesBetween returns all covering edges from→to (forward and inverse), in
// owner transformation order (accuracy ascending within each owner list).
func edgesBetween(from, to string, lon, lat float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) []graphEdge {
	ensureReverseIndex()
	from = strings.ToLower(from)
	to = strings.ToLower(to)

	var out []graphEdge

	glon := lon
	if fromDatum, err := loadDatum(from); err == nil {
		glon = lon + fromDatum.primeMeridianLongitude()
		for i, t := range fromDatum.Transformations {
			if t.Target == nil || t.Operation == nil {
				continue
			}

			if !timeSpecificAllowed(t.Operation, hasCoordEpoch, coordEpoch) {
				continue
			}

			if !strings.EqualFold(t.Target.Name, to) {
				continue
			}

			if !t.BoundingBox.Contains(glon, lat) {
				continue
			}

			k := edgeKey{owner: from, index: i, inverse: false}
			if excluded[k] {
				continue
			}

			out = append(out, graphEdge{key: k, to: to, accuracy: t.Accuracy, inverse: false})
		}
	}

	if toDatum, err := loadDatum(to); err == nil {
		// Inverse edge: coordinates remain in `from`'s frame (Greenwich AoU check).
		for i, t := range toDatum.Transformations {
			if t.Target == nil || t.Operation == nil {
				continue
			}

			if !timeSpecificAllowed(t.Operation, hasCoordEpoch, coordEpoch) {
				continue
			}

			if !strings.EqualFold(t.Target.Name, from) {
				continue
			}

			if !t.BoundingBox.Contains(glon, lat) {
				continue
			}

			k := edgeKey{owner: to, index: i, inverse: true}
			if excluded[k] {
				continue
			}

			out = append(out, graphEdge{key: k, to: to, accuracy: t.Accuracy, inverse: true})
		}
	}

	return out
}

type pathNode struct {
	name   string
	cost   float64
	hops   int
	viaHub bool // true if wgs84 appears as a non-terminal intermediate so far
	prev   int  // index in settled slice; -1 for start
	edge   graphEdge
}

type pqItem struct {
	idx  int
	cost float64
	hops int
	hub  bool
}

type pathPQ []pqItem

func (p pathPQ) Len() int {
	return len(p)
}

func (p pathPQ) Less(i, j int) bool {
	if p[i].cost != p[j].cost {
		return p[i].cost < p[j].cost
	}

	if p[i].hops != p[j].hops {
		return p[i].hops < p[j].hops
	}

	if p[i].hub != p[j].hub {
		return !p[i].hub && p[j].hub
	}

	return p[i].idx < p[j].idx
}
func (p pathPQ) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func (p *pathPQ) Push(x any) {
	*p = append(*p, x.(pqItem))
}

func (p *pathPQ) Pop() any {
	old := *p
	n := len(old)
	item := old[n-1]
	*p = old[:n-1]

	return item
}

func betterPath(aCost float64, aHops int, aHub bool, bCost float64, bHops int, bHub bool) bool {
	if aCost != bCost {
		return aCost < bCost
	}

	if aHops != bHops {
		return aHops < bHops
	}

	if aHub != bHub {
		return !aHub && bHub
	}

	return false
}

// findBestPath finds a minimum accumulated-accuracy path from→to covering (lon,lat).
func findBestPath(from, to Datum, lon, lat float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) ([]string, error) {
	ensureReverseIndex()

	start := datumName(from)
	goal := datumName(to)

	if start == "" || goal == "" {
		return nil, UnsupportedError{Err: errors.New("missing datum name")}
	}

	if start == goal {
		return []string{start}, nil
	}

	if excluded == nil {
		excluded = map[edgeKey]bool{}
	}

	// Datums with no published ops are treated as WGS84-equivalent (identity hops).
	searchFrom, searchTo := from, to
	prefix, suffix := "", ""

	if start != "wgs84" && isHubEquivalent(from) {
		prefix = start
		searchFrom = WGS84
	}

	if goal != "wgs84" && isHubEquivalent(to) {
		suffix = goal
		searchTo = WGS84
	}

	coreStart := datumName(searchFrom)
	coreGoal := datumName(searchTo)

	var core []string

	if coreStart == coreGoal {
		core = []string{coreStart}
	} else {
		var err error
		core, err = findBestPathCore(coreStart, coreGoal, lon, lat, excluded, hasCoordEpoch, coordEpoch)
		if err != nil {
			return nil, err
		}
	}

	path := make([]string, 0, len(core)+2)
	if prefix != "" {
		path = append(path, prefix)
	}

	path = append(path, core...)
	if suffix != "" {
		path = append(path, suffix)
	}

	return dedupeConsecutive(path), nil
}

func dedupeConsecutive(path []string) []string {
	if len(path) == 0 {
		return path
	}

	out := []string{path[0]}
	for _, n := range path[1:] {
		if n != out[len(out)-1] {
			out = append(out, n)
		}
	}

	return out
}

func findBestPathCore(start, goal string, lon, lat float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) ([]string, error) {
	nodes := []pathNode{{
		name: start,
		cost: 0,
		hops: 0,
		prev: -1,
	}}
	best := map[string]int{start: 0}

	pq := &pathPQ{{idx: 0, cost: 0, hops: 0, hub: false}}
	heap.Init(pq)

	for pq.Len() > 0 {
		item := heap.Pop(pq).(pqItem)
		cur := nodes[item.idx]
		if bestIdx, ok := best[cur.name]; ok && bestIdx != item.idx {
			continue
		}
		if cur.name == goal {
			return reconstructPath(nodes, item.idx), nil
		}
		if cur.hops >= maxPathHops {
			continue
		}

		for _, e := range outgoingEdges(cur.name, lon, lat, excluded, hasCoordEpoch, coordEpoch) {
			nextCost := cur.cost + e.accuracy
			nextHops := cur.hops + 1
			viaHub := cur.viaHub
			if cur.name == "wgs84" && start != "wgs84" {
				viaHub = true
			}

			if prevIdx, ok := best[e.to]; ok {
				prev := nodes[prevIdx]
				if !betterPath(nextCost, nextHops, viaHub, prev.cost, prev.hops, prev.viaHub) {
					continue
				}
			}

			idx := len(nodes)
			nodes = append(nodes, pathNode{
				name:   e.to,
				cost:   nextCost,
				hops:   nextHops,
				viaHub: viaHub,
				prev:   item.idx,
				edge:   e,
			})
			best[e.to] = idx
			heap.Push(pq, pqItem{idx: idx, cost: nextCost, hops: nextHops, hub: viaHub})
		}
	}

	return nil, UnsupportedError{Err: fmt.Errorf("no transformation path from %s to %s", start, goal)}
}

func reconstructPath(nodes []pathNode, idx int) []string {
	var rev []string
	for idx >= 0 {
		rev = append(rev, nodes[idx].name)
		idx = nodes[idx].prev
	}
	for i, j := 0, len(rev)-1; i < j; i, j = i+1, j-1 {
		rev[i], rev[j] = rev[j], rev[i]
	}
	return rev
}

func applyEdge(from Datum, e graphEdge, lon, lat, h float64) (float64, float64, float64, error) {
	owner, err := loadDatum(e.key.owner)
	if err != nil {
		return 0, 0, 0, err
	}
	if e.key.index < 0 || e.key.index >= len(owner.Transformations) {
		return 0, 0, 0, UnsupportedError{Err: errors.New("invalid transformation index")}
	}
	t := owner.Transformations[e.key.index]
	if e.inverse {
		return t.FromDatum(owner.Spheroid, lon, lat, h)
	}
	// Forward: owner must be `from`.
	_ = from
	return t.ToDatum(owner.Spheroid, lon, lat, h)
}

func applyHop(fromName, toName string, lon, lat, h float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) (float64, float64, float64, []edgeKey, error) {
	from, err := loadDatum(fromName)
	if err != nil {
		return 0, 0, 0, nil, err
	}
	toDatum, toErr := loadDatum(toName)

	cands := edgesBetween(fromName, toName, lon, lat, excluded, hasCoordEpoch, coordEpoch)
	if len(cands) == 0 {
		// Identity hop between WGS84 and a hub-equivalent datum (no published ops).
		if toErr == nil && identityHopAllowed(from, toDatum, fromName, toName) {
			return lon, lat, h, nil, nil
		}
		return 0, 0, 0, nil, UnsupportedError{
			Err: fmt.Errorf("no transformation from %s to %s", fromName, toName),
		}
	}

	var lastErr error
	tried := make([]edgeKey, 0, len(cands))
	for _, e := range cands {
		tried = append(tried, e.key)
		outLon, outLat, outH, err := applyEdge(from, e, lon, lat, h)
		if err != nil {
			lastErr = err
			continue
		}
		return outLon, outLat, outH, tried, nil
	}
	if lastErr == nil {
		lastErr = UnsupportedError{Err: errors.New("no valid transformation found")}
	}
	return 0, 0, 0, tried, lastErr
}

func identityHopAllowed(from, to Datum, fromName, toName string) bool {
	fromEq := fromName == "wgs84" || isHubEquivalent(from)
	toEq := toName == "wgs84" || isHubEquivalent(to)
	return fromEq && toEq
}

func applyPath(nodes []string, lon, lat, h float64, excluded map[edgeKey]bool, hasCoordEpoch bool, coordEpoch float64) (float64, float64, float64, []edgeKey, error) {
	if len(nodes) < 2 {
		return lon, lat, h, nil, nil
	}

	var failed []edgeKey

	for i := 0; i < len(nodes)-1; i++ {
		var tried []edgeKey
		var err error

		lon, lat, h, tried, err = applyHop(nodes[i], nodes[i+1], lon, lat, h, excluded, hasCoordEpoch, coordEpoch)

		failed = append(failed, tried...)
		if err != nil {
			return 0, 0, 0, failed, err
		}
	}

	return lon, lat, h, nil, nil
}
