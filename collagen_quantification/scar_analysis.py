# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import ndimage
import skimage as ski
import tifffile
from scipy import ndimage
import os, glob, time

# %% [markdown]
# > first, load in data

# %%
def load_pair(mask_path, tissue_path):
    """Return (collagen: bool array, tissue: bool array, px_um: float)."""

    mask = ski.io.imread(mask_path)
    tissue = ski.io.imread(tissue_path)

    mask_boolean = mask.astype(bool)
    tissue_boolean = tissue.astype(bool)

    with tifffile.TiffFile(mask_path) as mask_tif:
        xres = mask_tif.pages[0].tags['XResolution'].value
        px_um_mask = xres[1] / xres[0]
        unit_mask = (mask_tif.imagej_metadata or {}).get('unit')

    with tifffile.TiffFile(tissue_path) as tissue_tif:
        xres = tissue_tif.pages[0].tags['XResolution'].value
        px_um_tissue = xres[1] / xres[0]
        unit_tissue = (tissue_tif.imagej_metadata or {}).get('unit')

    valid_units = {'micron', 'microns', 'µm', 'um'}
    if unit_mask not in valid_units or unit_tissue not in valid_units:
        raise ValueError(f"{mask_path}: expected microns, got mask={unit_mask!r}, tissue={unit_tissue!r}")
    
    # Check that the mask and tissue images have the same shape and pixel size
    # then check that the pixel sizes are close enough to be considered equal
    # stop the function if they are not equal
    if mask_boolean.shape != tissue_boolean.shape:
        raise ValueError(f"Mask and tissue shapes do not match: {mask_boolean.shape} vs {tissue_boolean.shape}")

    if not np.isclose(px_um_mask, px_um_tissue):
        raise ValueError(f"Pixel sizes do not match: {px_um_mask} vs {px_um_tissue}")
    
    return mask_boolean, tissue_boolean, px_um_tissue # doesnt matter which pixel size we return since they are equal

# %%
image_path = "image_data/"

# %%
collagen_mask, tissue_mask, px_um = load_pair(
    os.path.join(image_path, "3135_mask.tif"),
    os.path.join(image_path, "3135_tissue.tif")
)

# %%
print(f"Collagen mask shape: {collagen_mask.shape}\nTissue shape: {tissue_mask.shape}")
print(f"Pixel size: {px_um} um")

# %% [markdown]
# > um conversion 0.325

# %%
# is all the collagen in the tissue? (i.e. is the collagen mask a subset of the tissue mask?)
np.sum((collagen_mask == True) & (tissue_mask == False))

# %%
# quick collagen area check
np.sum(collagen_mask) / np.sum(tissue_mask) * 100

# %% [markdown]
# > now calculate collagen density in a rolling window

# %%
def local_density(collagen, tissue, window_um, px_um):
    """Fraction of tissue that is collagen in a window around each pixel.
    Returns a float array, NaN outside tissue."""

    # convert window_um to pixels
    window_px = round(window_um / px_um)
    if window_px % 2 == 0:
        window_px += 1  # make it odd

    smooth_collagen_in_tissue = ndimage.uniform_filter(
        collagen & tissue,
        size = window_px,
        mode = 'constant',
        cval = 0,
        output = np.float32
        )
    
    smooth_tissue = ndimage.uniform_filter(
        tissue, 
        size = window_px, 
        mode = 'constant', 
        cval = 0,
        output = np.float32
        )

    density = np.divide(
        smooth_collagen_in_tissue, # collagen in tissue window
        smooth_tissue, # tissue window
        out=np.full_like(smooth_collagen_in_tissue, np.nan),
        where=smooth_tissue!=0
    )

    # report density only at tissue pixels (removes the halo around the tissue edge)
    density[~tissue.astype(bool)] = np.nan

    del smooth_collagen_in_tissue, smooth_tissue
    
    return density

# %%
res = local_density(collagen_mask, tissue_mask, window_um=100, px_um=px_um)

# %%
# downsample the result for visualization
small = ski.measure.block_reduce(res, block_size=10, func=np.nanmean)
plt.imshow(small, cmap='viridis')

# %%
print(f"Data type: {res[collagen_mask & tissue_mask].dtype}\nShape: {res[collagen_mask & tissue_mask].shape}\nMin: {np.nanmin(res[collagen_mask & tissue_mask])}, Max: {np.nanmax(res[collagen_mask & tissue_mask])}, Mean: {np.nanmean(res[collagen_mask & tissue_mask])}")

# %%
thr = np.quantile(res[collagen_mask & tissue_mask], 0.95)
print(f"95th percentile of collagen density in tissue: {thr}")

# %%
thr_mask = (res > thr) & collagen_mask & tissue_mask

# %%
small = ski.measure.block_reduce(thr_mask, block_size=10, func=np.max)
plt.imshow(small, cmap='Reds')

# %%
plt.hist(res[collagen_mask & tissue_mask], bins=100, range=(0, 1))
plt.show()

# %% [markdown]
# ---
# ## run full

# %%
image_path = "image_data/"

# %%
# find files and pairs by sample ID
collagen_mask_paths = sorted(glob.glob(os.path.join(image_path, "*_mask.tif")))
sample_IDs = [os.path.basename(f).split("_mask.tif")[0] for f in collagen_mask_paths]
tissue_mask_paths = [os.path.join(image_path, f"{s}_tissue.tif") for s in sample_IDs]

missing = [p for p in tissue_mask_paths if not os.path.exists(p)]
if missing:
    raise FileNotFoundError(f"Missing tissue masks: {missing}")

# settings here
windows_um = np.arange(50, 501, 50)          # sweep through 50, 100, ..., 500
bin_edges = np.linspace(0, 1, 1001)          # 1000 bins of width 0.001

# results
# store histograms for threshold cutoffs instead of raw values, which can be huge
n_img, n_win, n_bins = len(sample_IDs), len(windows_um), len(bin_edges) - 1
counts = np.zeros((n_img, n_win, n_bins), dtype=np.int64)
px_sizes = np.zeros(n_img)
collagen_pct = np.zeros(n_img)

# process one image at a time
for i, (mask_path, tissue_path, sample_id) in enumerate(
        zip(collagen_mask_paths, tissue_mask_paths, sample_IDs)):
    t0 = time.perf_counter()

    collagen, tissue, px_um = load_pair(mask_path, tissue_path)
    collagen_in_tissue = collagen & tissue

    px_sizes[i] = px_um
    collagen_pct[i] = 100 * collagen_in_tissue.sum() / tissue.sum()

    for j, w in enumerate(windows_um):
        density = local_density(collagen, tissue, w, px_um)
        counts[i, j], _ = np.histogram(density[collagen_in_tissue], bins=bin_edges)

    print(f"{sample_id}: collagen {collagen_pct[i]:.2f}%  ({time.perf_counter() - t0:.0f} s)")

# save results
np.savez("density_histograms.npz",
         sample_IDs=np.array(sample_IDs), windows_um=windows_um,
         bin_edges=bin_edges, counts=counts,
         px_sizes=px_sizes, collagen_pct=collagen_pct)

# %%
# load in results and calculate scar using a global percentile cutoff
results = np.load("density_histograms.npz")
print(results)

# %%
counts = results['counts'].copy()
windows_um = results['windows_um']
bin_edges = results['bin_edges']
print(f"histogram: {counts}\nDimensions: {counts.shape}")

# %% [markdown]
# > dimensions: number of samples, number of window sizes, number of observations (histogram bins)

# %%
# conver to 0 to 1 scale for percentile
norm_counts = counts / counts.sum(axis=-1, keepdims=True)
print(f"Normalized counts:\n{norm_counts}\nDimensions: {norm_counts.shape}")

# %%
# average normalized counts from all images within each distance window
averaged_norm_counts = norm_counts.mean(axis=0)
print(f"Average:\n{averaged_norm_counts}\nDimensions (should be n windows x n bins):\n{averaged_norm_counts.shape}")

# %%
# cumulative sum within each distance window
cumulative_dist = np.cumsum(averaged_norm_counts, axis=1)

# %%
# percentiles
percentiles = [0.85, 0.9, 0.95, 0.99]

# initialize array to store results
cutoff_table = np.full((len(windows_um), len(percentiles)), np.nan)
scar_bin_table = np.full((len(windows_um), len(percentiles)), -1, dtype=int)

# search along each row separately
for i, p in enumerate(percentiles):
    for j, w in enumerate(windows_um):
        first_scar_bin = np.searchsorted(cumulative_dist[j,], p) + 1 # find cutoff bin, add 1 so this is the upper edge
        scar_bin_table[j, i] = first_scar_bin
        cutoff_table[j, i] = bin_edges[first_scar_bin]

cutoff_df = pd.DataFrame(cutoff_table)
scar_bin_df = pd.DataFrame(scar_bin_table)

# name rows and columns
cutoff_df.columns = percentiles
cutoff_df.index = windows_um

scar_bin_df.columns = percentiles
scar_bin_df.index = windows_um

print(cutoff_df)
print(scar_bin_df)

# %%
n_images, n_windows, n_percentiles = norm_counts.shape[0], len(windows_um), len(percentiles)

scar_pct = np.full((n_images, n_windows, n_percentiles), np.nan)

for i, p in enumerate(percentiles):
    for j, w in enumerate(windows_um):
        scar_cutoff = scar_bin_table[j, i]
        # over all samples, sum at window j values from scar cutoff to maximum * 100 for scar%
        scar_pct[:, j, i] = norm_counts[:, j, scar_cutoff:].sum(axis=-1) * 100

# %%
sample_IDs = results['sample_IDs']         
windows_um = results['windows_um']          
collagen_pct = results['collagen_pct']     

idx = pd.MultiIndex.from_product(
    [sample_IDs, windows_um, percentiles],
    names=['sample_ID', 'window_um', 'percentile']
)

cutoff_full = np.broadcast_to(cutoff_table, scar_pct.shape)
collagen_full = np.broadcast_to(collagen_pct[:, None, None], scar_pct.shape)

scar_df = pd.DataFrame(
    {
        'cutoff': cutoff_full.ravel(),
        'scar_pct': scar_pct.ravel(),
        'collagen_pct': collagen_full.ravel(),
    },
    index=idx
).reset_index()

img, j, i = 2, 3, 1 
row = scar_df[(scar_df.sample_ID == sample_IDs[img]) &
              (scar_df.window_um == windows_um[j]) &
              (scar_df.percentile == percentiles[i])]
assert np.isclose(row.scar_pct.item(), scar_pct[img, j, i])
assert np.isclose(row.cutoff.item(), cutoff_table[j, i])
assert np.isclose(row.collagen_pct.item(), collagen_pct[img])

print(scar_df.shape)      
print(scar_df.head(8))

scar_df.to_csv("scar_results_long.csv", index=False)


