# Cluster Works

Interactive Streamlit isochrone fitting lab based on `Isochrone_Fitting.ipynb`.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Required data

The app requires three original binary/data files in the same directory as `app.py`:

- `isochrone_fitting.npz`: two synthetic observed cluster CMDs
- `parsec_ubv1.dat`: PARSEC UBV model grid (part 1)
- `parsec_ubv2.dat`: PARSEC UBV model grid (part 2)

**These large data files are not yet checked into this repository.** Copy them from the supplied `Cluster_Works.zip` package into the repository directory before running the app. The original `isochrones.py` is included to load the legacy serialized NumPy archive.

## Student activity

Students select one of two clusters, choose an age and metallicity from the original discrete model grid, and adjust the distance modulus to visually fit the B−V versus V CMD. No extinction correction or age/metallicity interpolation is applied.

The clusters' simulated ages are not exact grid points. The notebook is the source of the model grids and synthetic observations.
