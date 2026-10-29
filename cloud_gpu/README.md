# Cloud GPU Acceleration Guide: Google Colab & Kaggle

This directory contains GPU acceleration pipelines for running large protein language models (**ESM-2 650M** and **ProtBERT 420M**) on cloud accelerators (NVIDIA T4 / V100 / A100), bypassing local CPU memory limits and slow download speeds.

---

## Option 1: Running in Google Colab (Recommended)

### Using the Colab Extension in Antigravity / VS Code:
1. Click the **Colab** icon in the sidebar (highlighted in your editor).
2. Click **Sign in with Google** if not already signed in.
3. Open `cloud_gpu/colab_hsv_plm_benchmark.ipynb`.
4. In the top right of the notebook, select **Google Colab** as the kernel and choose **GPU (T4)**.
5. Click **Run All**.
6. Weights (1.68 GB) will download in ~10 seconds over Google Cloud's 10 Gbps connection, and full 929-mutation masked marginal scoring will finish in ~60 seconds!

### Using Web Browser Google Colab:
1. Go to [colab.research.google.com](https://colab.research.google.com).
2. Click **Upload** -> Select `cloud_gpu/colab_hsv_plm_benchmark.ipynb`.
3. In the menu: **Runtime** -> **Change runtime type** -> Select **T4 GPU** -> Save.
4. Drag and drop `primary_position_folds.tsv` and `hsv_plm_gpu_runner.py` into the Files pane (`/content/`).
5. Run all cells.

---

## Option 2: Running via Kaggle CLI

Because Kaggle CLI credentials (`~/.kaggle/access_token`) are already verified on this machine:

1. Push and launch the GPU kernel:
   ```powershell
   kaggle kernels push -p d:\antiviral-resistance-hsv--1\cloud_gpu
   ```
2. Check execution status:
   ```powershell
   kaggle kernels status <username>/hsv-plm-gpu-benchmark
   ```
3. Retrieve output results:
   ```powershell
   kaggle kernels output <username>/hsv-plm-gpu-benchmark -p d:\antiviral-resistance-hsv--1\logs
   ```
