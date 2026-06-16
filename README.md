# Voronoi Analysis

Computes a per-atom Voronoi tessellation and reports the topological coordination number, the atomic volume, and the cavity radius (the largest empty sphere centred on the atom that fits inside its cell).

## Install

```bash
vpm install @voltlabs/voronoi-analysis
```

## CLI

```bash
voronoi-analysis <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--cutoff <float>` | no | `0` | Neighbor search cutoff radius. `0` requests an automatic estimate from the mean atomic volume. |
| `--edgeThreshold <float>` | no | `0` | Minimum edge length for a face to contribute to the coordination count. |
| `--faceThreshold <float>` | no | `0` | Minimum face area for a face to contribute to the coordination count. |
| `--useRadii` | no | off | Build a radical (power) Voronoi tessellation using the `Radius` column. |
| `--onlySelected` | no | off | Only tessellate atoms whose `Selection` column is non-zero. |
| `--threads <int>` | no | auto | Maximum worker threads (OneTBB). |
| `--help` | no | — | Print CLI help. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_voronoi.parquet` | Voronoi Analysis | — *(listing-only)* |
| `{output_base}_atoms.parquet` | Voronoi Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/voronoi-analysis
