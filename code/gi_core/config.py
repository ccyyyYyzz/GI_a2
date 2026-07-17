"""Frozen global constants — ROUND59 spec §0/§1. Do not modify after preregistration."""

SEED0 = 20260717

M = 20000          # frames per (family, seed) run
MW = 200000        # witness-only pattern draws (no bucket)
SEEDS_A = [0, 1, 2]
SEEDS_B = [0, 1, 2, 3, 4]
SEEDS_EXT = [0, 1, 2]

PHOTONS_B = [1e3, 1e4]
PHOTONS_A = [1e4]

N_A = 256    # 16x16, Phase A
N_B = 4096   # 64x64, Phase B

# family ids (stable stream ids for RNG derivation — frozen)
FAMILY_IDS = {
    "GAUSS": 1,
    "GAM1": 2,
    "GAM2": 3,
    "GAM3": 4,
    "GAM4": 5,
    "GAM8": 6,
    "CORR-LOGN": 7,
    "MIX-LOGN": 8,
}

LINK_IDS = {
    "LIN": 1,
    "DT30": 2,
    "FGAIN": 3,
    "SAT30": 4,
    "SAT50": 5,
    "GAMMA07": 6,
    "LOG": 7,
}

# stream tags
STREAM_PATTERN = 0
STREAM_NOISE = 1
STREAM_WITNESS = 2
STREAM_PROBES = 3
STREAM_ESTIMATOR = 4

# GAUSS / CORR-LOGN spatial kernel: correlation length 2 px, marginal sigma
CORR_LEN_PX = 2.0
GAUSS_SIGMA = 0.25
LOGN_SIGMA2 = 0.25          # sigma_ln^2
MIX_CORR_LENS = (2.0, 3.0)  # two states, same marginal
KERNEL_JITTER_REL = 1e-8    # relative diagonal jitter; part of the model definition
                            # (used consistently in sampling, score, and WHITEN-OR)

READOUT_COEF = 2.0          # sigma_r = READOUT_COEF / s
DT30_TAU_COEF = 3.0 / 7.0   # tau = (3/7)/s  ->  tau*lambda_mean = 3/7, mean compression 0.7
FGAIN_CV = 0.2

# metrics
IMAGES = ["cat", "deer", "dog", "horse", "airplane", "car", "ship", "truck"]

# gates (preregistered)
A1_MEDIAN_FACTOR = 1.5
A1_P90_FACTOR = 2.0
A3_PSNR_DB = 0.5
A3_MIN_POS = 6
B1_PSNR_DB = 0.50
B1_SSIM = 0.010
B1_LPIPS = 0.010
B1_MIN_POS = 6
BOOT_N = 10000
BOOT_LB_Q = 0.10
