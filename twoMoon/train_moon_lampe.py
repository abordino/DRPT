# train_moon_lampe_multi.py
# LAMPE NRE/BNRE on the moon simulator for budgets {2^6, 2^9, 2^12, 2^15}.
# Saves TorchScript models: artifacts/moon_rhat_{nre|bnre}_B{B}.pt

import os
import torch
import torch.nn as nn
import torch.optim as optim
import zuko
from lampe.inference import NRE, NRELoss
from lampe.inference.bnre import BNRELoss
from lampe.utils import GDStep
from lampe.data import JointDataset

# ---------------- 0) Device & seeds ----------------
device = "cuda" if torch.cuda.is_available() else "cpu"
torch.manual_seed(42)

# ---------------- 1) Prior & simulator -------------
LOWER = -torch.ones(2)
UPPER =  torch.ones(2)
prior = zuko.distributions.BoxUniform(LOWER, UPPER)

def simulator(theta: torch.Tensor) -> torch.Tensor:
    n = theta.shape[:-1]
    alpha = (torch.rand(*n, device=theta.device) - 0.5) * torch.pi      # ~ Unif(-pi/2, pi/2)
    rr    = 0.01 * torch.randn(*n, device=theta.device) + 0.1           # ~ N(0.1, 0.01^2)
    base  = torch.stack([rr * torch.cos(alpha) + 0.25, rr * torch.sin(alpha)], dim=-1)
    shift = torch.stack([
        -torch.abs(theta[..., 0] + theta[..., 1]) / (2.0 ** 0.5),
        (-theta[..., 0] + theta[..., 1]) / (2.0 ** 0.5),
    ], dim=-1)
    return base + shift

# ---------------- 2) Helpers -----------------------
def build_dataset(B: int, batch_size: int = 256) -> JointDataset:
    # sample once (fixed simulation budget)
    theta = prior.sample((B,))                       # (B, 2) on CPU
    x     = simulator(theta)                         # (B, 2)
    return JointDataset(theta, x, batch_size=batch_size, shuffle=True)

def train(estimator: NRE, loss_module, dataset: JointDataset,
          epochs: int = 400, lr: float = 1e-3, device: str = device):
    est = estimator.to(device)
    opt = optim.AdamW(est.parameters(), lr=lr)
    step = GDStep(opt, clip=1.0)
    est.train()
    for _ in range(epochs):
        for theta_b, x_b in dataset:                 # dataset yields (batch,2) tensors on CPU
            theta_b = theta_b.to(device)
            x_b     = x_b.to(device)
            step(loss_module(theta_b, x_b))
    est.eval()
    return est

class RatioWrapper(torch.nn.Module):
    def __init__(self, est):
        super().__init__()
        self.est = est       # returns log r
    def forward(self, theta: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return torch.exp(self.est(theta, x))  # r_hat = exp(log r)

def export_torchscript(estimator: NRE, out_path: str):
    # move to CPU for portable TorchScript
    wrapper = RatioWrapper(estimator.to("cpu")).eval()
    example_th, example_x = torch.zeros(1, 2), torch.zeros(1, 2)
    ts = torch.jit.trace(wrapper, (example_th, example_x))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    ts.save(out_path)

# ---------------- 3) Main loop over budgets --------
if __name__ == "__main__":
    budgets = [4]  # 64, 512, 4096, 32768
    for B in budgets:
        print(f"\n=== Training with budget B = {B} ===")
        trainset = build_dataset(B, batch_size=256)

        # Fresh estimators per budget
        est_nre  = NRE(2, 2, hidden_features=[128]*5, activation=nn.ELU)
        est_bnre = NRE(2, 2, hidden_features=[128]*5, activation=nn.ELU)

        # Losses
        loss_nre  = NRELoss(est_nre)
        loss_bnre = BNRELoss(est_bnre, lmbda=100.0)

        # Train
        est_nre  = train(est_nre,  loss_nre,  trainset, epochs=400, lr=1e-3)
        est_bnre = train(est_bnre, loss_bnre, trainset, epochs=400, lr=1e-3)

        # Export
        export_torchscript(est_nre,  f"artifacts/moon_rhat_nre_B{B}.pt")
        export_torchscript(est_bnre, f"artifacts/moon_rhat_bnre_B{B}.pt")

        print(f"Saved:")
        print(f"  artifacts/moon_rhat_nre_B{B}.pt")
        print(f"  artifacts/moon_rhat_bnre_B{B}.pt")
