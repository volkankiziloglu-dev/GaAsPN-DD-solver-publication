# Third-party input data (not redistributed)

Place the third-party files listed in `docs/data_sources.md` in this folder, using the file
names (and, for the digitized Klitzke curves, the sub-folder `experimental/klitzke2025/`) given
there. Alternatively set the environment variable `GAASPN_EXTERNAL_DATA` to a folder with the
same layout. `src/resolveDataFile.m` searches `src/data/`, this folder and
`$GAASPN_EXTERNAL_DATA`, in this order.

Everything in this folder except this README is ignored by git.
