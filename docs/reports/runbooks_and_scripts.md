# IAM Runbooks & Server/Agent Startup Scripts

> **Note:** Only the first 10 search results are shown. [View more results in GitHub](https://github.com/lewispoul/IAM/search?q=start+run+watch+serve+agent+noxctl+ctl+runner+background+daemon.sh+server+background+watch+port+env+python+flask+daemon+run+start&type=code)

---

## 1. Scripts That Start Servers, Agents, or Watchers

### A. Typical Startup Scripts
- **Names searched:** `start_all.sh`, `run_server.sh`, `watcher.sh`, `noxctl`, `runner.sh`, `agent.sh`, etc.
- **Common locations:** Project root, `IAM_GUI/`, `tools/`, `scripts/`, or `bin/`.

### B. Example Usage Patterns

#### a. Starting Flask Web Server
```bash
export FLASK_APP=IAM_GUI/backend.py
export FLASK_ENV=production
export IAM_PORT=5000  # default
python3 IAM_GUI/backend.py
# or
flask run --host=0.0.0.0 --port=${IAM_PORT:-5000}
```

#### b. Starting Agent/Watcher
```bash
python3 IAM_Agent.py &
# Watches for new jobs in notebooks/results/ or IAM_Knowledge/
export IAM_AGENT_MODE=batch
export IAM_AGENT_WATCH_DIR=notebooks/results
```

#### c. Starting All Services
```bash
./start_all.sh
# May invoke: backend.py, agent.py, watcher scripts, UI server, etc.
```

---

## 2. Environment Variables

| Variable           | Purpose                                           | Example Value      |
|--------------------|--------------------------------------------------|--------------------|
| FLASK_APP          | Main backend entrypoint                          | IAM_GUI/backend.py |
| FLASK_ENV          | Flask environment (dev/prod)                     | production         |
| IAM_PORT           | Port for web server                              | 5000, 5006         |
| IAM_AGENT_MODE     | Agent operation mode                             | batch, watcher     |
| IAM_AGENT_WATCH_DIR| Directory for file watching                      | notebooks/results  |

---

## 3. Ports Used

- **Flask backend:** `5000` (default), `5006` (WSL or alternate)
- **Agent/Watcher:** background, no direct port unless serving API
- **UI frontend:** usually same as backend, e.g., `http://localhost:5000`

---

## 4. Expected Folder Layout

```
IAM/
├── IAM_GUI/
│   ├── backend.py          # Flask server main
│   ├── static/             # JS, CSS, frontend assets
│   ├── templates/          # HTML UI templates
│   ├── iam_fixed.js        # Robust frontend logic
│   └── performance_utils.py
├── IAM_Agent.py            # Batch/async agent
├── IAM_Knowledge/          # Results (JSON, CSV)
├── notebooks/
│   └── results/            # Per-job files, watched by agent
├── tools/
│   ├── xtb/                # XTB binaries, test files
│   ├── psi4/               # Psi4 test files
├── scripts/
│   ├── start_all.sh        # (if present) starts all services
│   └── watcher.sh          # (if present) watches for changes
```

---

## 5. Summary Table

| Script/File      | Usage                      | Env Vars         | Port        | Folders Involved          |
|------------------|---------------------------|------------------|-------------|--------------------------|
| backend.py       | Flask web server          | FLASK_APP, FLASK_ENV | 5000    | IAM_GUI/, static/, templates/ |
| IAM_Agent.py     | Background agent/watcher  | IAM_AGENT_MODE   | N/A         | notebooks/results/, IAM_Knowledge/ |
| start_all.sh     | Start all services        | Various          | Multiple    | Project root              |

**Results may be incomplete due to search limits.**
