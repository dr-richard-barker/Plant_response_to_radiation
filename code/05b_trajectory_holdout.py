#!/usr/bin/env python3
"""
05b_trajectory_holdout.py
Phase 2.2b — Within-study time-trajectory holdout evaluation.

The most scientifically meaningful test of the GP-AE: for each study with a
time series, hold out the LATEST timepoint(s), train on the earlier ones, and
test whether the model predicts the late-timepoint response. This directly
evaluates the user's goal of "predicting gene expression trajectories under
varying radiation stress scenarios."

We also generate trajectory samples along a dense (dose, time, LET) grid using
the trained decoder, producing the Phase 2.2 deliverable: a continuous
pseudo-temporal expression curve.
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

import importlib.util
spec = importlib.util.spec_from_file_location("gp_ae", "/mnt/results/zenodo_bundle/code/05_gp_autoencoder.py")
gp_ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(gp_ae)

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
OUT_DIR = Path("/mnt/results/zenodo_bundle/results/trajectories")
CKPT_DIR = Path("/mnt/shared-workspace/shared/checkpoints")
DEVICE = gp_ae.DEVICE


def _safe_spearman(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if a.size < 3 or b.size < 3: return float("nan")
    if np.nanstd(a) < 1e-12 or np.nanstd(b) < 1e-12: return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r, _ = spearmanr(a, b)
    return float(r) if not np.isnan(r) else float("nan")


def trajectory_holdout(X, C, obs, epochs=100):
    """For each study with a time series, hold out the latest timepoint and
    test whether the model predicts it from earlier timepoints.

    We report two metrics:
      - R^2 on centered expression (absolute reconstruction; expected to be
        modest because extrapolation to an unseen timepoint is hard).
      - Direction-of-change correlation: Spearman rho between the predicted
        response vector (x_hat - previous_timepoint_mean) and the actual
        response vector (x_test - previous_timepoint_mean). This measures
        whether the model gets the *direction* of the late response right,
        which is the scientifically meaningful question.
    """
    results = []
    studies = obs["StudyID"].values
    times = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce").values
    X_np = X.cpu().numpy()

    print("\n=== Within-study time-trajectory holdout ===")
    print("(hold out latest timepoint, predict from earlier ones)")
    for sid in np.unique(studies):
        if not sid: continue
        mask = studies == sid
        study_times = times[mask]
        unique_times = np.unique(study_times[~np.isnan(study_times)])
        if len(unique_times) < 2: continue
        max_time = unique_times.max()
        # Second-to-last timepoint as the reference for direction-of-change
        if len(unique_times) >= 2:
            ref_time = np.sort(unique_times)[-2]
        else:
            ref_time = unique_times[0]

        global_train = mask & np.array([(t != max_time and not np.isnan(t)) for t in times])
        global_test = mask & (times == max_time)
        global_ref = mask & (times == ref_time)
        if global_train.sum() < 5 or global_test.sum() < 1: continue

        model, _ = gp_ae.train_model(X[global_train], C[global_train], X.shape[1],
                                      epochs=epochs, verbose=False)
        res = gp_ae.evaluate(model, X[global_test], C[global_test], obs[global_test])
        n_train, n_test = int(global_train.sum()), int(global_test.sum())

        # Direction-of-change: predicted vs actual response relative to ref timepoint
        with torch.no_grad():
            x_hat_test, _, _, _, _ = model(X[global_test].to(DEVICE), C[global_test].to(DEVICE))
            x_hat_test = x_hat_test.cpu().numpy()
        x_ref = X_np[global_ref].mean(axis=0) if global_ref.sum() > 0 else X_np[global_train].mean(axis=0)
        x_test = X_np[global_test]
        actual_delta = (x_test - x_ref).mean(axis=0)  # mean response at late time
        pred_delta = (x_hat_test - x_ref).mean(axis=0)
        dir_rho = _safe_spearman(actual_delta, pred_delta)
        # Also cosine similarity (direction alignment)
        cos_sim = float(np.dot(actual_delta, pred_delta) /
                        (np.linalg.norm(actual_delta)*np.linalg.norm(pred_delta) + 1e-8))

        print(f"  {sid}: hold out t={max_time}h (n_test={n_test}, n_train={n_train})  "
              f"R^2={res['R2']:.3f}  dir_rho={dir_rho:.2f}  cos_sim={cos_sim:.2f}")
        results.append({"study": str(sid), "held_time_h": float(max_time),
                        "ref_time_h": float(ref_time),
                        "n_train": n_train, "n_test": n_test,
                        "R2": res["R2"], "rho_time": res["rho_time"],
                        "direction_rho": float(dir_rho if not np.isnan(dir_rho) else 0.0),
                        "cosine_sim": cos_sim})
    return results


def sample_trajectories(model, obs, gene_names, cov_mean, cov_std, gene_mean, gene_std):
    """Sample expression along a dense (dose, time, LET) grid using the decoder.
    Produces the continuous pseudo-temporal expression curve deliverable."""
    model.eval()
    # Build a grid: dose x time x LET, varying one covariate at a time around
    # representative values, so we get interpretable 1D trajectories.
    # Representative center: median dose, median time, median LET.
    log_dose_center = 0.0  # already standardized
    log_time_center = 0.0
    log_let_center = 0.0

    # Dose trajectory: vary log_dose from -2 to +2 std, fix time & LET at center
    n_pts = 25
    dose_grid = np.linspace(-2, 2, n_pts)
    time_grid = np.linspace(-2, 2, n_pts)
    let_grid = np.linspace(-2, 2, n_pts)

    trajectories = {}
    with torch.no_grad():
        # Use a fixed latent (the cohort mean latent) to isolate covariate effect
        # Encode all samples, take mean latent
        X_all, C_all, _, _, _, _, _, _ = gp_ae.prepare_data()
        mu_all, _ = model.encoder(X_all.to(DEVICE))
        z_mean = mu_all.mean(dim=0)  # shape [latent_dim]
        z_mean = z_mean.unsqueeze(0).expand(n_pts, -1).to(DEVICE)  # [n_pts, latent_dim]

        for name, grid in [("dose", dose_grid), ("time", time_grid), ("let", let_grid)]:
            covs = np.zeros((n_pts, 3))
            idx = {"dose": 0, "time": 1, "let": 2}[name]
            covs[:, idx] = grid
            # Other two covariates at center (0)
            covs_t = torch.tensor(covs, dtype=torch.float32).to(DEVICE)
            x_hat = model.decoder(z_mean, covs_t).cpu().numpy()
            # Un-standardize: x_hat is in standardized centered-expression space
            x_unnorm = x_hat * gene_std + gene_mean  # back to centered log1p
            trajectories[name] = x_unnorm

    # Save as long-format CSV: Covariate, Value, Gene, Expression
    rows = []
    # Convert grid values back to original units for interpretability
    for name, grid in [("dose", dose_grid), ("time", time_grid), ("let", let_grid)]:
        idx = {"dose": 0, "time": 1, "let": 2}[name]
        orig_vals = grid * cov_std[idx] + cov_mean[idx]
        orig_vals = np.exp(orig_vals)  # back from log space
        x = trajectories[name]
        for i, val in enumerate(orig_vals):
            for g, gene in enumerate(gene_names):
                rows.append({"covariate": name, "value": float(val),
                             "gene": gene, "expression": float(x[i, g])})
    traj_df = pd.DataFrame(rows)
    traj_path = OUT_DIR / "decoder_trajectories.csv"
    traj_df.to_csv(traj_path, index=False)
    print(f"Decoder trajectories -> {traj_path}  ({len(traj_df)} rows)")
    return traj_df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Device: {DEVICE}")
    print("Preparing data...")
    X, C, obs, gene_names, cov_mean, cov_std, gene_mean, gene_std = gp_ae.prepare_data()
    print(f"  Expression: {X.shape}")

    # Trajectory holdout evaluation
    traj_results = trajectory_holdout(X, C, obs, epochs=100)
    with open(OUT_DIR / "trajectory_holdout.json", "w") as f:
        json.dump(traj_results, f, indent=2, default=float)
    if traj_results:
        mean_r2 = np.mean([r["R2"] for r in traj_results])
        print(f"\nTrajectory holdout mean R^2: {mean_r2:.3f}")

    # Load trained model and sample decoder trajectories
    print("\nLoading trained model for trajectory sampling...")
    ckpt = torch.load(CKPT_DIR / "gp_ae_checkpoint.pt", map_location=DEVICE, weights_only=False)
    model = gp_ae.GPAutoencoder(X.shape[1], gp_ae.LATENT_DIM, n_covariates=3).to(DEVICE)
    model.load_state_dict(ckpt["model_state"])
    model.to(DEVICE)
    traj_df = sample_trajectories(model, obs, gene_names, cov_mean, cov_std,
                                   gene_mean, gene_std)

    # Also save latent colored by covariates for visualization
    print("\nExtracting latent embeddings for visualization...")
    model.eval()
    with torch.no_grad():
        _, _, _, z, _ = model(X.to(DEVICE), C.to(DEVICE))
        z = z.cpu().numpy()
    z_df = pd.DataFrame(z, index=obs["SampleID"], columns=[f"z{i}" for i in range(gp_ae.LATENT_DIM)])
    z_df["StudyID"] = obs["StudyID"].values
    z_df["RadiationQuality"] = obs["RadiationQuality"].values
    z_df["AbsorbedDose_Gy"] = pd.to_numeric(obs["AbsorbedDose_Gy"], errors="coerce").values
    z_df["TimePostExposure_h"] = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce").values
    z_df["LET_keV_um"] = pd.to_numeric(obs["LET_keV_um"], errors="coerce").values
    z_df.to_csv(OUT_DIR / "latent_embeddings_annotated.csv")
    print(f"Annotated latent -> {OUT_DIR / 'latent_embeddings_annotated.csv'}")


if __name__ == "__main__":
    main()
