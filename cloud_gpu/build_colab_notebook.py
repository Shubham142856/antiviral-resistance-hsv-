import json

nb = {
    "cells": [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# HSV-1 Antiviral Resistance: Cloud GPU Protein Language Model Benchmark\n",
                "\n",
                "**Project:** HSV Antiviral Resistance Mutation Prediction (Paper 2)  \n",
                "**Models:** ESM-2 (35M, 650M) and ProtBERT (420M)  \n",
                "**Standard:** Dual-Split Benchmark (Position-Grouped CV vs. Contiguous-Block CV)  \n",
                "**Stopping Rule:** $\\ge 0.65$ AUROC with 95% Position-Cluster Bootstrap CI $> 0.50$ under Contiguous-Block holdout on $\\ge 2/3$ primary pairs (`UL23-ACV`, `UL30-ACV`, `UL30-FOS`).\n",
                "\n",
                "> **How to run in Google Colab:**\n",
                "> 1. Runtime -> Change runtime type -> **T4 GPU** (or A100 GPU).\n",
                "> 2. Click **Run all** (Ctrl+F9).\n",
                "> 3. Results will be printed, visualized, and saved to `plm_benchmark_results.json`."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# 1. Environment Setup & Hardware Verification\n",
                "!nvidia-smi\n",
                "!pip install -q transformers torch pandas numpy scikit-learn matplotlib seaborn"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# 2. Dataset Preparation (Download or Load Clean Single-Site Cohort)\n",
                "import os\n",
                "import pandas as pd\n",
                "\n",
                "data_file = 'primary_position_folds.tsv'\n",
                "if not os.path.exists(data_file):\n",
                "    print('Downloading primary_position_folds.tsv...')\n",
                "    # Fallback to direct raw github link or local mount\n",
                "    !curl -s -O https://raw.githubusercontent.com/ojcharles/herpesdrg-db/main/inst/extdata/primary_position_folds.tsv || true\n",
                "\n",
                "if os.path.exists(data_file):\n",
                "    df = pd.read_csv(data_file, sep='\\t')\n",
                "    print(f'Successfully loaded {len(df)} primary cohort rows across {df[\"gene\"].nunique()} genes.')\n",
                "else:\n",
                "    print('Please upload primary_position_folds.tsv to the runtime /content folder.')"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# 3. Run Benchmark Runner Script\n",
                "!python hsv_plm_gpu_runner.py"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# 4. Display Results & Generate Visualizations\n",
                "import json\n",
                "import matplotlib.pyplot as plt\n",
                "import seaborn as sns\n",
                "\n",
                "if os.path.exists('plm_benchmark_results.json'):\n",
                "    with open('plm_benchmark_results.json', 'r') as f:\n",
                "        results = json.load(f)\n",
                "    \n",
                "    print(json.dumps(results, indent=2))\n",
                "    \n",
                "    # Visual comparison bar chart\n",
                "    pairs = ['UL23-Aciclovir', 'UL30-Aciclovir', 'UL30-Foscarnet']\n",
                "    fig, ax = plt.subplots(1, 2, figsize=(14, 5))\n",
                "    \n",
                "    for m_key, data in results.items():\n",
                "        if 'metrics' not in data: continue\n",
                "        pos_means = [data['metrics'][p]['pos_grouped_mean'] for p in pairs]\n",
                "        blk_means = [data['metrics'][p]['contig_block_mean'] for p in pairs]\n",
                "        ax[0].plot(pairs, pos_means, marker='o', label=m_key)\n",
                "        ax[1].plot(pairs, blk_means, marker='s', label=m_key)\n",
                "        \n",
                "    ax[0].set_title('Position-Grouped CV (Weak Split)')\n",
                "    ax[0].set_ylim(0.4, 1.0)\n",
                "    ax[0].axhline(0.65, color='gray', linestyle='--', label='Threshold 0.65')\n",
                "    ax[0].grid(True, alpha=0.3)\n",
                "    ax[0].legend()\n",
                "    \n",
                "    ax[1].set_title('Contiguous-Block CV (Strict Structural Split)')\n",
                "    ax[1].set_ylim(0.4, 1.0)\n",
                "    ax[1].axhline(0.65, color='red', linestyle='--', label='Stopping Threshold 0.65')\n",
                "    ax[1].grid(True, alpha=0.3)\n",
                "    ax[1].legend()\n",
                "    \n",
                "    plt.tight_layout()\n",
                "    plt.savefig('plm_benchmark_summary.png', dpi=300)\n",
                "    plt.show()\n",
                "    print('Saved plot to plm_benchmark_summary.png')\n",
                "else:\n",
                "    print('plm_benchmark_results.json not found.')"
            ]
        }
    ],
    "metadata": {
        "accelerator": "GPU",
        "colab": {
            "gpuType": "T4",
            "provenance": []
        },
        "kernelspec": {
            "display_name": "Python 3",
            "name": "python3"
        },
        "language_info": {
            "name": "python"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 0
}

with open(r'd:\antiviral-resistance-hsv--1\cloud_gpu\colab_hsv_plm_benchmark.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=2)

print('Generated colab_hsv_plm_benchmark.ipynb!')
