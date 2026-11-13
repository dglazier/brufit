# Reweight Simulations to Match Data (Importance Sampling)

This example demonstrates how to use the `HS::FIT::Weights` class in **brufit** to reweight a simulated Monte Carlo (MC) distribution so that it matches a target data distribution via importance sampling.

---

## Quick Start

### 1. Generate Mock Data and Simulation
Run `CreateData.C` with ROOT to generate toy datasets:

```bash
root -l -q CreateData.C
```

This creates two files with $10^5$ events each, containing an event ID branch (`fgID`) and a photon energy branch (`Eg`):
* **`data.root`** (`TTree` `"data"`): `Eg` sampled from a broader Gaussian, $\mathcal{N}(\mu=9, \sigma=10)$.
* **`sim.root`** (`TTree` `"sim"`): `Eg` sampled from a narrower Gaussian, $\mathcal{N}(\mu=9, \sigma=7)$.

### 2. Compute Importance-Sampling Weights
Run `CreateWeights.C` using `brufit`:

```bash
brufit CreateWeights.C
```

This script performs the following steps:
1. Loads the `"data"` tree from `data.root` and the `"sim"` tree from `sim.root`.
2. Configures a `Weights` object named `"MCWeights"` with species name `"LikeData"`, event ID branch `"fgID"`, and output file `"Weights.root"`.
3. Defines a binning template histogram (`TH1F`) and calls `ImportanceSampling()` to compute the Data/MC ratio as a function of `Eg` and assign a weight to each simulated event.
4. Sorts the weight index by event ID (`SortWeights()`), plots the unweighted MC, target data, and reweighted MC (`Draw1DWithWeights()`) on the same canvas for visual verification, and saves the weights to `Weights.root`.

---

## Reweighting to Background-Subtracted Data (sWeights)

If your target data has existing weights (for example, `Signal` sWeights from a prior `sPlot` fit) or your MC already has prior weights, pass those `Weights` objects into `ImportanceSampling()`:

```cpp
// Load existing Data sWeights
Weights* dataWeights = new Weights();
dataWeights->LoadSavedDisc("path/to/Tweights.root", "HSsWeights");

// Compute MC weights matching the "Signal"-weighted data distribution
impWeights->ImportanceSampling(treeMC, treeData, weightHist, "Eg",
                               nullptr, "",             // Optional prior MC weights & species
                               dataWeights, "Signal");  // Target Data weights & species
```

---

## Using the MC Weights in a BruFit Fit

Once `Weights.root` is saved, you can pass these weights into the MC integration of any `BruEventsPDF` (including `BruComponentsPDF`) by appending `"WEIGHTS@<Species>,<File>,<ObjectName>"` to the PDF factory or parser string:

* **Via `Setup::FactoryPDF` or `Setup::ParserPDF`:**
  Append `"WEIGHTS@LikeData,Weights.root,MCWeights"` to the end of the PDF configuration string.
* **Important:** Ensure `Fitter.SetUp().SetIDBranchName("fgID");` in your fit script matches the ID branch name used when generating `Weights.root`.