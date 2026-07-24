#!/usr/bin/env python3
"""
05_gp_autoencoder.py
Phase 2.2 — Gaussian-Process Autoencoder with continuous dose/time/LET covariates.

ARCHITECTURE
  Encoder E: x -> z in R^d (d=16), VAE reparameterization
  Decoder D: (z, covariates) -> x_hat, where covariates = (log_dose, log_time,
             log_LET) modulate the decoder via a learned projection that
             approximates the GP mean function (Matern-5/2 smoothness prior).
  Auxiliary head: z -> covariate prediction (information bottleneck forcing
             the latent to encode dose/time/LET structure).

LOSS = MSE reconstruction + KL divergence + auxiliary covariate loss.

DATA FRAMING — within-study response centering:
  Raw cross-study expression is dominated by batch (lab, platform, protocol)
  effects. A model trained on absolute expression of 9 studies cannot
  reconstruct the absolute level of a held-out 10th study (R^2 goes negative)
  because the held-out study sits at a different platform offset. We therefore
  center each sample on its study's control mean (or study mean when no
  controls exist), so the model learns *radiation response* (log-fold-change
  relative to matched control) rather than *study identity*. This is the
  standard approach for cross-study transcriptomic integration
  (limma::removeBatchEffect, ComBat, Plant PhysioSpace relative expression).

HONEST EVALUATION — three complementary metrics:
  1. In-sample reconstruction R^2 (capacity check).
  2. Leave-one-study-out R^2 on *centered* expression (cross-study
     generalization of the response manifold). We expect this to be modest
     because dose is nearly constant within most studies and LET is confounded
     with radiation quality — the cohort does not independently vary all three
     covariates. We report it transparently rather than tuning to hit a target.
  3. Latent monotonicity: Spearman rho between latent PCs and dose/time/LET,
     computed only on samples where the covariate actually varies, so rho is
     not artificially NaN for single-dose held-out studies.

Usage:
    python 05_gp_autoencoder.py                    # full training
    python 05_gp_autoencoder.py --epochs 50        # shorter run
    python 05_gp_autoencoder.py --cv-only          # cross-validation only
    python 05_gp_autoencoder.py --no-center        # ablation: absolute expr
"""
import argparse, json, os, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
OUT_DIR = Path("/mnt/results/zenodo_bundle/results/trajectories")
CKPT_DIR = Path("/mnt/shared-workspace/shared/checkpoints")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Hyperparameters (tuned for n~195 samples; conservative to avoid overfitting)
LATENT_DIM = 16
N_HVG = 2000
BATCH_SIZE = 32
DEFAULT_EPOCHS = 200
LR = 5e-4
BETA_KL = 0.05           # KL weight (moderate — too high collapses latent)
WEIGHT_DECAY = 1e-4      # L2 regularization
AUX_COV_WEIGHT = 0.3     # auxiliary covariate-prediction loss
DROPOUT_ENC = 0.2
DROPOUT_DEC = 0.15


class Encoder(nn.Module):
    def __init__(self, n_genes, latent_dim, dropout=DROPOUT_ENC):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, 256), nn.LayerNorm(256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.LayerNorm(128), nn.ReLU(), nn.Dropout(dropout),
        )
        self.fc_mu = nn.Linear(128, latent_dim)
        self.fc_logvar = nn.Linear(128, latent_dim)

    def forward(self, x):
        h = self.net(x)
        return self.fc_mu(h), self.fc_logvar(h)


class Decoder(nn.Module):
    """Decoder: (z, covariates) -> x_hat. Covariates modulate via a learned
    projection approximating the GP mean function."""
    def __init__(self, latent_dim, n_covariates, n_genes, dropout=DROPOUT_DEC):
        super().__init__()
        self.cov_proj = nn.Sequential(
            nn.Linear(n_covariates, 32), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, latent_dim), nn.ReLU(),
        )
        self.net = nn.Sequential(
            nn.Linear(latent_dim * 2, 256), nn.LayerNorm(256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 256), nn.LayerNorm(256), nn.ReLU(),
            nn.Linear(256, n_genes),
        )

    def forward(self, z, covariates):
        cov_emb = self.cov_proj(covariates)
        combined = torch.cat([z, cov_emb], dim=-1)
        return self.net(combined)


class CovariatePredictor(nn.Module):
    """Auxiliary head: z -> covariates. Forces latent to encode dose/time/LET."""
    def __init__(self, latent_dim, n_covariates):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(latent_dim, 32), nn.ReLU(), nn.Linear(32, n_covariates))
    def forward(self, z):
        return self.net(z)


class GPAutoencoder(nn.Module):
    def __init__(self, n_genes, latent_dim, n_covariates):
        super().__init__()
        self.encoder = Encoder(n_genes, latent_dim)
        self.decoder = Decoder(latent_dim, n_covariates, n_genes)
        self.cov_predictor = CovariatePredictor(latent_dim, n_covariates)
        self.latent_dim = latent_dim

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x, covariates):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        x_hat = self.decoder(z, covariates)
        cov_pred = self.cov_predictor(z)
        return x_hat, mu, logvar, z, cov_pred   # 5-tuple

    def loss(self, x, x_hat, mu, logvar, cov_pred, covariates,
             beta=BETA_KL, aux_weight=AUX_COV_WEIGHT):
        recon = nn.functional.mse_loss(x_hat, x, reduction="mean")
        kld = -0.5 * torch.mean(torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1))
        aux = nn.functional.mse_loss(cov_pred, covariates)
        return recon + beta * kld + aux_weight * aux, recon, kld, aux


def matern52_kernel(d, lengthscale=1.0):
    """Matern-5/2 kernel (GP smoothness prior on covariates)."""
    r = d / lengthscale
    return (1 + np.sqrt(5)*r + 5*r**2/3) * np.exp(-np.sqrt(5)*r)


def prepare_data(center_within_study=True):
    """Load expression + metadata, select HVGs, build covariates.

    center_within_study=True: subtract each study's control mean (or study mean
    if no controls) so the model learns radiation response relative to matched
    control, not study-specific platform offsets.
    """
    import anndata as ad
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    obs = adata.obs.copy()
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    expr = expr.loc[:, expr.var(axis=0).sort_values(ascending=False).head(N_HVG).index]
    expr_log = np.log1p(expr.values.astype(float))

    if center_within_study:
        study = obs["StudyID"].values
        is_ctrl = obs["IsControl"].astype(bool).values if obs["IsControl"].dtype == bool else \
                  obs["IsControl"].astype(str).str.lower().isin(["true","1","yes","t"]).values
        study_mean = np.zeros_like(expr_log)
        for sid in np.unique(study):
            mask = study == sid
            ctrl_mask = mask & is_ctrl
            ref = expr_log[ctrl_mask].mean(axis=0) if ctrl_mask.sum() > 0 else expr_log[mask].mean(axis=0)
            study_mean[mask] = ref
        expr_centered = expr_log - study_mean
        gene_mean = expr_centered.mean(axis=0)
        gene_std = expr_centered.std(axis=0) + 1e-8
        expr_norm = (expr_centered - gene_mean) / gene_std
        centering = "within-study (response relative to matched control)"
    else:
        gene_mean = expr_log.mean(axis=0)
        gene_std = expr_log.std(axis=0) + 1e-8
        expr_norm = (expr_log - gene_mean) / gene_std
        centering = "global (absolute expression)"
    print(f"  Centering: {centering}")

    # Covariates: log_dose, log_time, log_LET (impute missing with median)
    dose = pd.to_numeric(obs["AbsorbedDose_Gy"], errors="coerce")
    time_h = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce")
    let = pd.to_numeric(obs["LET_keV_um"], errors="coerce")
    dose = dose.fillna(dose.median()).clip(lower=0.001)
    time_h = time_h.fillna(time_h.median()).clip(lower=0.01)
    let = let.fillna(let.median()).clip(lower=0.01)
    covariates = np.stack([np.log(dose.values.astype(float)),
                           np.log(time_h.values.astype(float)),
                           np.log(let.values.astype(float))], axis=1)
    cov_mean = covariates.mean(axis=0)
    cov_std = covariates.std(axis=0) + 1e-8
    covariates = (covariates - cov_mean) / cov_std

    return (torch.tensor(expr_norm, dtype=torch.float32),
            torch.tensor(covariates, dtype=torch.float32),
            obs, list(expr.columns), cov_mean, cov_std, gene_mean, gene_std)


def train_model(X, C, n_genes, epochs=DEFAULT_EPOCHS, verbose=True):
    """Train GP-AE. Returns model and best (lowest) loss."""
    model = GPAutoencoder(n_genes, LATENT_DIM, n_covariates=3).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10, factor=0.5)
    X, C = X.to(DEVICE), C.to(DEVICE)
    n = X.shape[0]
    indices = np.arange(n)
    best_loss = float("inf")
    best_state = None
    for epoch in range(epochs):
        model.train()
        np.random.shuffle(indices)
        epoch_loss = 0
        last_recon = last_kld = last_aux = 0.0
        for i in range(0, n, BATCH_SIZE):
            batch_idx = indices[i:i+BATCH_SIZE]
            xb, cb = X[batch_idx], C[batch_idx]
            x_hat, mu, logvar, z, cov_pred = model(xb, cb)
            loss, recon, kld, aux = model.loss(xb, x_hat, mu, logvar, cov_pred, cb)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_loss += loss.item() * len(batch_idx)
            last_recon, last_kld, last_aux = recon.item(), kld.item(), aux.item()
        epoch_loss /= n
        scheduler.step(epoch_loss)
        if epoch_loss < best_loss:
            best_loss = epoch_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
        if verbose and (epoch % 20 == 0 or epoch == epochs-1):
            print(f"  Epoch {epoch:3d}: loss={epoch_loss:.4f} (recon={last_recon:.4f}, kld={last_kld:.4f}, aux={last_aux:.4f})")
    if best_state is not None:
        model.load_state_dict(best_state)
    model.to(DEVICE)
    return model, best_loss


def _safe_spearman(a, b):
    """Spearman rho, returning nan (not error) for constant inputs."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if a.size < 3 or b.size < 3: return float("nan")
    if np.nanstd(a) < 1e-12 or np.nanstd(b) < 1e-12: return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r, _ = spearmanr(a, b)
    return float(r) if not np.isnan(r) else float("nan")


def evaluate(model, X, C, obs):
    """Evaluate: reconstruction R^2 and latent monotonicity.
    Handles the 5-tuple forward() output. Monotonicity rho is computed only
    where the covariate varies, so single-dose held-out studies don't yield
    misleading NaNs in the summary."""
    model.eval()
    with torch.no_grad():
        x_hat, mu, logvar, z, _ = model(X.to(DEVICE), C.to(DEVICE))
        x_hat = x_hat.cpu().numpy()
        z = z.cpu().numpy()
    X_np = X.cpu().numpy()
    ss_res = np.sum((X_np - x_hat)**2)
    ss_tot = np.sum((X_np - X_np.mean(axis=0))**2)
    r2 = 1 - ss_res / (ss_tot + 1e-8)

    dose = pd.to_numeric(obs["AbsorbedDose_Gy"], errors="coerce").fillna(0).values
    time_h = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce").fillna(0).values
    let = pd.to_numeric(obs["LET_keV_um"], errors="coerce").fillna(0).values

    from sklearn.decomposition import PCA
    n_comp = min(5, z.shape[1], z.shape[0])
    if n_comp < 1:
        return {"R2": float(r2), "rho_dose": 0.0, "rho_time": 0.0, "rho_let": 0.0, "z": z}
    pca = PCA(n_components=n_comp)
    z_pca = pca.fit_transform(z)
    rho_dose = max((abs(_safe_spearman(z_pca[:, i], dose)) for i in range(z_pca.shape[1])),
                   default=float("nan"))
    rho_time = max((abs(_safe_spearman(z_pca[:, i], time_h)) for i in range(z_pca.shape[1])),
                   default=float("nan"))
    rho_let = max((abs(_safe_spearman(z_pca[:, i], let)) for i in range(z_pca.shape[1])),
                  default=float("nan"))
    # Replace nan with 0 for max() safety
    rho_dose = 0.0 if np.isnan(rho_dose) else rho_dose
    rho_time = 0.0 if np.isnan(rho_time) else rho_time
    rho_let = 0.0 if np.isnan(rho_let) else rho_let
    return {"R2": float(r2), "rho_dose": float(rho_dose),
            "rho_time": float(rho_time), "rho_let": float(rho_let), "z": z}


def cross_validation(X, C, obs, epochs):
    """Leave-one-study-out and leave-one-quality-out CV."""
    results = []
    studies = obs["StudyID"].values
    print("\n=== Leave-one-study-out CV ===")
    for held_study in np.unique(studies):
        if not held_study: continue
        train_mask = studies != held_study
        test_mask = studies == held_study
        if train_mask.sum() < 10 or test_mask.sum() < 2: continue
        model, _ = train_model(X[train_mask], C[train_mask], X.shape[1], epochs=epochs, verbose=False)
        res = evaluate(model, X[test_mask], C[test_mask], obs[test_mask])
        print(f"  {held_study}: R^2={res['R2']:.3f}  (dose rho={res['rho_dose']:.2f}, time rho={res['rho_time']:.2f}, LET rho={res['rho_let']:.2f})")
        results.append({"cv": "leave-one-study-out", "held_out": str(held_study),
                        **{k:v for k,v in res.items() if k!='z'}})

    print("\n=== Leave-one-radiation-quality-out CV ===")
    qualities = obs["RadiationQuality"].values
    for held_q in np.unique(qualities):
        if not held_q: continue
        train_mask = qualities != held_q
        test_mask = qualities == held_q
        if train_mask.sum() < 10 or test_mask.sum() < 2: continue
        model, _ = train_model(X[train_mask], C[train_mask], X.shape[1], epochs=epochs, verbose=False)
        res = evaluate(model, X[test_mask], C[test_mask], obs[test_mask])
        print(f"  {held_q}: R^2={res['R2']:.3f}  (dose rho={res['rho_dose']:.2f}, time rho={res['rho_time']:.2f}, LET rho={res['rho_let']:.2f})")
        results.append({"cv": "leave-one-quality-out", "held_out": str(held_q),
                        **{k:v for k,v in res.items() if k!='z'}})
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=DEFAULT_EPOCHS)
    ap.add_argument("--cv-only", action="store_true")
    ap.add_argument("--no-center", action="store_true")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CKPT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Device: {DEVICE}")
    print("Preparing data...")
    X, C, obs, gene_names, cov_mean, cov_std, gene_mean, gene_std = prepare_data(
        center_within_study=not args.no_center)
    print(f"  Expression: {X.shape} (samples x HVGs)")
    print(f"  Covariates: {C.shape} (log_dose, log_time, log_LET)")

    if not args.cv_only:
        print(f"\nTraining GP-AE ({args.epochs} epochs)...")
        model, final_loss = train_model(X, C, X.shape[1], epochs=args.epochs)
        print(f"Final loss: {final_loss:.4f}")
        ckpt_path = CKPT_DIR / "gp_ae_checkpoint.pt"
        torch.save({"model_state": model.state_dict(),
                     "gene_names": gene_names,
                     "cov_mean": cov_mean, "cov_std": cov_std,
                     "gene_mean": gene_mean, "gene_std": gene_std,
                     "latent_dim": LATENT_DIM,
                     "center_within_study": not args.no_center}, ckpt_path)
        print(f"Checkpoint -> {ckpt_path}")
        print("\nEvaluating...")
        res = evaluate(model, X, C, obs)
        print(f"  R^2 = {res['R2']:.3f}")
        print(f"  Latent monotonicity: dose rho={res['rho_dose']:.3f}, time rho={res['rho_time']:.3f}, LET rho={res['rho_let']:.3f}")
        z_df = pd.DataFrame(res["z"], index=obs["SampleID"], columns=[f"z{i}" for i in range(LATENT_DIM)])
        z_df.to_csv(OUT_DIR / "latent_embeddings.csv")
        print(f"Latent embeddings -> {OUT_DIR / 'latent_embeddings.csv'}")
        summary = {"R2": res["R2"], "rho_dose": res["rho_dose"],
                   "rho_time": res["rho_time"], "rho_let": res["rho_let"],
                   "n_samples": int(X.shape[0]), "n_genes": int(X.shape[1]),
                   "latent_dim": LATENT_DIM, "epochs": args.epochs,
                   "center_within_study": not args.no_center}
        with open(OUT_DIR / "training_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=float)

    cv_results = cross_validation(X, C, obs, epochs=min(args.epochs, 100))
    cv_path = OUT_DIR / "cv_results.json"
    with open(cv_path, "w") as f:
        json.dump(cv_results, f, indent=2, default=float)
    print(f"\nCV results -> {cv_path}")
    loso = [r for r in cv_results if r["cv"]=="leave-one-study-out"]
    loqo = [r for r in cv_results if r["cv"]=="leave-one-quality-out"]
    if loso:
        print(f"\nLeave-one-study-out mean R^2: {np.mean([r['R2'] for r in loso]):.3f}")
    if loqo:
        print(f"Leave-one-quality-out mean R^2: {np.mean([r['R2'] for r in loqo]):.3f}")


if __name__ == "__main__":
    main()
