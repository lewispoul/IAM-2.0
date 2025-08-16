# IAM Unified Environment Bill of Materials

> **Note:** Only the first 10 search results are shown. For more, [see full code search results in GitHub](https://github.com/lewispoul/IAM/search?q=requirements.txt+Dockerfile+Makefile+.sh+.yaml+.yml+conda+environment+pip+install+cuda+gpu+nvidia+torch+rdkit+psi4+xtb&type=code).

---

## 1. Python Requirements (`requirements.txt` or pip)

- **Pinned core packages:**  
  - `rdkit==2023.03.1` *(chemistry toolkit, CPU only)*  
  - `pytest==7.4.0` *(testing)*

- **GPU options:**  
  - `torch==2.2.2` *(for ML, requires CUDA if GPU used)*  
  - `cupy==12.3.0` *(for GPU numpy/scipy drop-in, optional)*

---

## 2. Conda Environment (`environment.yml`)

```yaml
name: iam-env
channels:
  - conda-forge
dependencies:
  - python=3.10
  - rdkit=2023.03.1
  - psi4=1.8.0
  - numpy=1.23.5
  - scipy=1.10.1
  - pandas=2.0.3
  - matplotlib=3.7.2
  - flask=2.3.2
  - joblib=1.3.2
  - requests=2.31.0
  - xtb-python=0.3.0
```

**Notes:**  
- Add `cudatoolkit=11.7` for GPU support in conda.
- Psi4 and RDKit are available via conda-forge.
- XTB must be installed via system package or its own installer.

---

## 3. Dockerfile

```dockerfile
FROM nvidia/cuda:12.2.0-base-ubuntu22.04

RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/
WORKDIR /app
RUN pip install --no-cache-dir -r requirements.txt

# For GPU ML
RUN pip install torch==2.2.2 cupy==12.3.0

EXPOSE 5000
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "backend:app"]
```

**Notes:**  
- Use `nvidia/cuda` base image for GPU.  
- Psi4 and XTB may require manual installation/scripts if not apt-gettable.

---

## 4. Makefile (if present)

```makefile
env:
	conda env create -f environment.yml

run:
	gunicorn --bind 0.0.0.0:5000 backend:app

test:
	pytest tests/
```

---

## 5. Shell Scripts (`install.sh`, etc.)

```bash
#!/bin/bash
conda create -n iam-env python=3.10 rdkit=2023.03.1 psi4=1.8.0 numpy=1.23.5 scipy=1.10.1 pandas=2.0.3 matplotlib=3.7.2 flask=2.3.2 joblib=1.3.2 requests=2.31.0 xtb-python=0.3.0 -c conda-forge
conda activate iam-env
pip install torch==2.2.2 cupy==12.3.0
```

---

## 6. System Packages

- **XTB:** Must be installed from [https://github.com/grimme-lab/xtb](https://github.com/grimme-lab/xtb).  
  - Build instructions: `make`, `make install`
- **Psi4:** `conda install psi4 -c conda-forge`
- **GPU:** Use `nvidia-driver`, `cudatoolkit`, and compatible CUDA version.

---

## 7. Optional ML/Benchmark Dependencies

- `scikit-learn==1.3.0`
- `xgboost==2.0.0`
- `seaborn==0.12.2`

---

## 8. Summary Table

| Package        | Version    | CPU/GPU | Source           |
|----------------|------------|---------|------------------|
| python         | 3.10       | CPU     | conda/pip/Docker |
| rdkit          | 2023.03.1  | CPU     | conda/pip        |
| psi4           | 1.8.0      | CPU/GPU | conda/pip        |
| xtb            | latest     | CPU     | system/manual    |
| torch          | 2.2.2      | GPU     | pip/conda        |
| cupy           | 12.3.0     | GPU     | pip              |
| numpy          | 1.23.5     | CPU/GPU | conda/pip        |
| flask          | 2.3.2      | CPU     | conda/pip        |
| pandas         | 2.0.3      | CPU     | conda/pip        |

**Results may be incomplete due to search limits.**
