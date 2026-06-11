# VoronoiAnalysis

`VoronoiAnalysis` computes a per-atom Voronoi tessellation and reports the
topological coordination number, the atomic volume and the cavity radius
(the largest empty sphere centred on the atom that fits inside its cell).

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- VoronoiAnalysis
```

## CLI

Usage:

```bash
voronoi-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--cutoff <float>` | No | Neighbor search cutoff radius. `0` requests an automatic estimate from the mean atomic volume. | `0` |
| `--edgeThreshold <float>` | No | Minimum edge length for a face to contribute to the coordination count. | `0` |
| `--faceThreshold <float>` | No | Minimum face area for a face to contribute to the coordination count. | `0` |
| `--useRadii` | No | Build a radical (power) Voronoi tessellation using the `Radius` column. | off |
| `--onlySelected` | No | Only tessellate atoms whose `Selection` column is non-zero. | off |
| `--threads <int>` | No | Maximum worker threads (OneTBB). | auto |
| `--help` | No | Print CLI help. | |

### Outputs

Two Parquet files are written alongside `output_base`:

- `{output_base}_voronoi.parquet` - summary table (single-row JSON payload)
  with the `main_listing` statistics (`total_atoms`, `mean_atomic_volume`,
  `mean_coordination`, `effective_cutoff`, thresholds, `use_radii`).
- `{output_base}_atoms.parquet` - the canonical per-atom table consumed by
  coloring/filtering, the GLB export and listings. Fixed columns
  (`atom_index`, `id`, `x/y/z`, `bucket`, `structure_id`, `structure_name`,
  `cluster_id`) plus the Voronoi metrics: `coordination`, `atomic_volume`
  and `cavity_radius`. Atoms are grouped by `bucket = Coordination_<k>` with
  `structure_id = k`, driving the `AtomisticExporter` GLB.

## Algorithm

This plugin ships a self-contained Voronoi implementation based on
iterative half-space clipping of a convex polytope:

1. For each atom, collect all neighbors inside `cutoff` (either
   user-supplied or auto-estimated as `4 * cbrt(3 V_atom / 4pi)`).
2. Seed the cell as an axis-aligned cube of half-side `cutoff` centred on
   the atom.
3. For every neighbor (closest first) clip the cell against the
   perpendicular bisector of the neighbor vector (or the radical plane
   when `--useRadii` is requested).
4. After all clips, compute the cell volume by fan-triangulating each
   face and summing signed tetrahedra from the origin, and the cavity
   radius as the minimum origin-to-face-plane distance.
5. Count surviving faces that pass the area and edge-length filters to
   get the topological coordination number.

> **Note on topology accuracy.** This implementation is fast, portable
> (zero external geometry dependencies) and accurate for the vast
> majority of MD environments, including crystalline, liquid and
> amorphous systems. It does **not** implement the full numerical
> degeneracy handling of Rycroft's voro++ library. If you need
> publication-grade topology statistics for near-degenerate geometries,
> use OVITO's voro++-backed Voronoi modifier. For typical production MD
> trajectories the results from this plugin are in close agreement with
> voro++.
