import h5py

filename = "input/flat-shell/v7.3/matv73.mat"
with h5py.File(filename, "r") as f:
    K = f["K"]
    print("K keys:", list(K.keys()))
    for k in K.keys():
        print(k, K[k].shape, K[k].dtype)
