import numpy as np, torch, sys
print("PY", sys.version.split()[0], "np", np.__version__, "torch", torch.__version__,
      "cuda", torch.cuda.is_available(),
      "dev", (torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"))
import fog_e6 as fe
x = fe.make_scene(32); U = fe.dct_basis(48, 32)
print("SANITY_OK scene", x.shape, "U", U.shape, "orthoerr",
      float(np.abs(U.T @ U - np.eye(48)).max()))
