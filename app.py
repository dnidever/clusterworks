"""Cluster Works: interactive PARSEC isochrone fitting for stellar astrophysics."""
from pathlib import Path
import __main__
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from astropy.table import Table, vstack
from isochrones import Isochrone

ROOT = Path(__file__).resolve().parent
AGES = [0.1, 0.5, 1.0, 2.5, 5.0, 7.5, 10.0]
METALS = [-2.0, -1.5, -1.0, -0.5, 0.0]

st.set_page_config(page_title="Cluster Works", page_icon="✨", layout="wide")

@st.cache_resource(show_spinner="Loading cluster observations and PARSEC models…")
def load_data():
    # The supplied NPZ was saved from a script where Isochrone was in __main__.
    # Register the original class there to support this legacy pickle format.
    __main__.Isochrone = Isochrone
    with np.load(ROOT / "isochrone_fitting.npz", allow_pickle=True) as archive:
        obs = {key: archive[key].item() for key in ("obs1", "obs2")}
    # Rebuild the model grid from the supplied original PARSEC tables, avoiding
    # reliance on the embedded pickled model and preserving notebook selections.
    tables = [Table.read(ROOT / name, format="ascii") for name in
              ("parsec_ubv1.dat", "parsec_ubv2.dat")]
    tab = vstack(tables)
    age = 10.0 ** np.asarray(tab["logAge"], float) / 1e9
    mh = np.asarray(tab["MH"], float)
    keep = np.asarray(tab["label"]) < 8
    models = {}
    for a in AGES:
        for m in METALS:
            mask = keep & (np.abs(age - a) < 0.01) & (np.abs(mh - m) < 0.01)
            subset = tab[mask]
            models[(a, m)] = (np.asarray(subset["Bmag"], float),
                              np.asarray(subset["Vmag"], float))
    return obs, models

try:
    observations, models = load_data()
except Exception as exc:
    st.error(f"Unable to load the supplied files: {exc}")
    st.stop()

st.title("✨ Cluster Works")
st.caption("Explore stellar cluster ages, metallicities, and distances using PARSEC isochrones")
st.markdown("Fit a theoretical isochrone to the observed **B−V versus V** color–magnitude diagram. "
            "Change the model age, metallicity, and distance modulus until the main-sequence "
            "turnoff and evolved stars line up with the observations.")

with st.sidebar:
    st.header("Cluster and model")
    cluster = st.radio("Observed cluster", ["Cluster 1", "Cluster 2"])
    age = st.select_slider("Age (Gyr)", options=AGES, value=0.5)
    metal = st.select_slider("Metallicity [M/H]", options=METALS, value=-1.5)
    mu = st.slider("Distance modulus (mag)", min_value=5.0, max_value=22.0,
                   value=11.0 if cluster == "Cluster 1" else 16.0, step=0.1)
    show_iso = st.checkbox("Show theoretical isochrone", value=True)
    show_points = st.checkbox("Show model points", value=False)
    max_stars = st.slider("Observed stars displayed", 500, 50000, 12000, 500)
    st.divider()
    st.markdown("**Distance modulus**")
    st.latex(r"m-M=5\log_{10}(d/10\,\mathrm{pc})")
    st.metric("Implied distance", f"{10 ** ((mu + 5) / 5) / 1000:,.2f} kpc")

obs = observations["obs1" if cluster == "Cluster 1" else "obs2"]
b = np.asarray(obs["bmag"], dtype=float)
v = np.asarray(obs["vmag"], dtype=float)
valid = np.flatnonzero(np.isfinite(b) & np.isfinite(v))
# Stable random subsample so model changes never shuffle the observations.
rng = np.random.default_rng(372)
if len(valid) > max_stars:
    valid = np.sort(rng.choice(valid, size=max_stars, replace=False))

left, right = st.columns([1.65, 1], gap="large")
with left:
    st.subheader(f"{cluster}: Color–magnitude diagram")
    fig, ax = plt.subplots(figsize=(5.4, 5.7), constrained_layout=True)
    ax.scatter(b[valid] - v[valid], v[valid], s=2, alpha=0.28,
               color="#526477", rasterized=True, label="Observed stars")
    bi, vi = models[(age, metal)]
    good = np.isfinite(bi) & np.isfinite(vi)
    if show_iso and np.any(good):
        ax.plot((bi - vi)[good], vi[good] + mu, color="#e4572e", lw=2,
                label=f"PARSEC: {age:g} Gyr, [M/H]={metal:+.1f}")
        if show_points:
            ax.scatter((bi - vi)[good], vi[good] + mu, s=5,
                       color="#e4572e", alpha=0.65)
    elif show_iso:
        st.warning("No PARSEC points were found for this age and metallicity.")
    ax.set_xlabel("B − V (mag)")
    ax.set_ylabel("Apparent V magnitude")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim((10, 25) if cluster == "Cluster 1" else (15, 26))
    ax.invert_yaxis() if not ax.yaxis_inverted() else None
    ax.grid(alpha=0.16)
    ax.legend(loc="best", fontsize=9, markerscale=3)
    st.pyplot(fig, clear_figure=True, use_container_width=False)
with right:
    st.subheader("Your current fit")
    st.metric("Age", f"{age:g} Gyr")
    st.metric("Metallicity", f"{metal:+.1f} dex")
    st.metric("Distance modulus", f"{mu:.1f} mag")
    st.metric("Distance", f"{10 ** ((mu + 5) / 5) / 1000:.2f} kpc")
    st.markdown("### What should you look for?")
    st.markdown("**Age:** Focus on the main-sequence turnoff and the shape of the evolved-star sequence.\n\n"
                "**Metallicity:** Watch the color and morphology of the main sequence and giant branch.\n\n"
                "**Distance:** Shift the model vertically until the absolute-magnitude sequence "
                "matches the observed apparent magnitudes.")
    st.info("The model grid is discrete in age and metallicity, as in the notebook. "
            "Only the distance modulus changes continuously. Extinction is not fitted.")

with st.expander("Student questions / lab activity"):
    st.markdown("1. Fit Cluster 1. Record your age, [M/H], distance modulus, and distance.\n"
                "2. Repeat for Cluster 2. Which cluster is older and which is more metal-poor?\n"
                "3. What CMD features most strongly constrain age?\n"
                "4. Try changing only the metallicity. Can a change in age compensate?\n"
                "5. How does changing distance differ from changing age or metallicity?\n"
                "6. Explain two sources of uncertainty or degeneracy in your fits.")
with st.expander("About the data"):
    st.write("The two observed clusters are synthetic datasets supplied with the "
             "Isochrone Fitting notebook. The model curves come from the supplied "
             "PARSEC UBV tables, with pulsating AGB entries (label ≥ 8) excluded. "
             "The plotted models use absolute B and V magnitudes shifted by distance modulus. "
             "The notebook's cluster-generating script specifies approximate ages of "
             "0.4 and 6.5 Gyr, which do not necessarily coincide with the discrete "
             "model ages available to students.")
