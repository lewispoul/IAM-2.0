<!--
Copilot: This document is the master Implementation & Test Plan.
Follow the steps PHASE 0→5 and TEST PLAN A0→D1 in order.
After completing each task:
- Ensure Acceptance Criteria (AC) are satisfied.
- Run the test instructions (see TEST PLAN).
- Only proceed if tests pass.
-->
IAM 2.0 – Implementation & Test Plan (FULL)
0) Objectifs & principes

But: une plate‑forme modulaire pour simuler, prédire et visualiser (moteurs XTB/Psi4, prédicteurs VoD/Pcj, spectres, orbitales, stabilité…), en conservant les acquis d’IAM 1.0 tout en modernisant l’architecture (FastAPI + schéma d’E/S unifié + tests systématiques).

Principes

Un schéma de résultats unifié (clé/valeur + unités + artefacts) s’applique à tous les moteurs et pipelines.

Chaque feature possède : dépendances, endpoints, AC, tests rapides (curl/httpie), et marqueurs pytest.

Les tâches sont ordonnées et feature‑flagged pour livrer incrémentalement.

1) Standards d’E/S et persistance (héritage IAM 1.0 → 2.0)
1.1 Unified Results Schema — [Spec-ID: IO.1] 🔑

Portée

Modèles Pydantic (v2) pour toutes réponses:
schema_version, molecule_id, source, inputs, metrics, units, artifacts, created_at

Tous metrics avec unités explicites (eV, D, cm⁻¹, GPa, m/s…).

artifacts: chemins sûrs vers logs, cubes, spectres (CSV/JSON), PNG, ZIP.

AC

 Toute réponse passe la validation Pydantic.

 Version de schéma figée et émise.

 Rejet si une unité manque.

Smoke test (conceptuel)

pytest -k test_schema_conformance -q

1.2 Results Writer — [Spec-ID: IO.2] 🔑

Portée

Dossiers immuables par run: IAM_Knowledge/Results/<calc_id>/

Fichiers: result.json, provenance.json (+ artefacts).

Helper env/paths: IAM_RESULTS_BASE, IAM_EXPORTS_BASE.

AC

 calc_id unique ; write‑once.

 provenance.json contient versions moteurs + temps mur.

Pré‑requis

backend/utils/env.py, backend/utils/persistence.py (faits si non présents).

2) Pipelines de base (héritage IAM 1.0 + extensions)
2.1 SMILES→XYZ — [Spec-ID: PIPE.SMILES] 🔑

Portée

RDKit ETKDG + UFF/MMFF (embedd + opt légère).

Validateur XYZ (nombre d’atomes, formats).

AC

 SMILES invalides → 400 + envelope d’erreur standard.

 XYZ conforme → utilisé par moteurs.

Test

curl -X POST /api/convert/to-xyz -d '{"smiles":"CCO"}'

2.2 OPT→PROP (cache) — [Spec-ID: PIPE.OPT]

Portée

Hachage de job (input + options) ; réutilisation géométrie convergée.

Abstraction “properties from optimized geometry”.

AC

 Cache hit détectable ; gain temps mesuré.

 Résultats reproductibles (tolérances documentées).

Test

Deux runs identiques → 2nd plus rapide, mêmes metrics ± tolérance.

2.3 Orchestrateur de prédiction — [Spec-ID: PIPE.PRED]

Portée

Normalisation unités → appel prédicteurs (KJ, Keshavarz, ML).

Chaînage vers panneau “Performance”.

AC

 Entrées normalisées (densité, ΔHf, etc.).

 Sorties avec incertitudes si dispo.

Test

curl -X POST /api/v1/predict_vod -d '{...}' → unités correctes.

3) Moteurs de chimie (héritage + nouvelles intégrations)
3.1 XTB Wrapper (SP/OPT/FREQ) — [Spec-ID: ENG.XTB] 🔑

Portée

SP/OPT/FREQ réels (flag IAM_ENABLE_XTB), fallback stub sinon.

Extraction: énergie (hartree/eV), dipôle (vecteur + total), gap HOMO‑LUMO, fréquences (cm⁻¹), cubes (optionnel).

Timeout + kill d’arbre de processus ; logs propres.

AC

 Mode réel commutable par flag ; stub par défaut OK tests.

 Unités cohérentes ; artefacts écrits ; erreurs mappées (504 timeout, 500 moteur).

Tests

pytest -k xtb_stub (sans flag); IAM_ENABLE_XTB=1 pytest -k xtb_real.

3.2 Psi4 Wrapper (SP/OPT/FREQ+Spectres) — [Spec-ID: ENG.PSI4]

Portée

Intégration locale (flag IAM_ENABLE_PSI4).

freq → spectres IR (intensités) (+ Raman si activé).

Sorties dans schéma unifié + spectres JSON/CSV.

AC

 SP/OPT/FREQ OK ; IR JSON conforme (cm⁻¹, intensité).

 Erreurs mappées (422 inputs, 500 moteur).

Test

IAM_ENABLE_PSI4=1 curl -X POST /api/v1/run_psi4 -d '{...,"calc":"freq"}'

3.3 Gaussian HPC Stub — [Spec-ID: ENG.GAUS]

Portée

Génération d’input .gjf + empaquetage ; pas d’exécution locale.

Connecteur “queue” distant (SSH / API), à implémenter plus tard.

AC

 Artefacts: .gjf, manifest, zip prêt à soumettre.

 Zéro exécution locale (sécurité/ressources).

Test

curl -X POST /api/v1/gaussian/prepare -d '{...}' → zip.

4) Prédicteurs énergétiques (héritage IAM 1.0 + 2.0)
4.1 Kamlet–Jacobs — [Spec-ID: PRED.KJ] 🔑

Portée

VoD (m/s), Pcj (GPa) à partir de ρ et Q.

Constants et hypothèses documentées.

AC

 Unités SI toujours en sortie ; validation de domaine.

 Cas limites (ρ manquante, Q manquant) → 422.

Test

curl -X POST /api/v1/predict_vod -d '{ "formula":"CH3NO2", "density":1.14, "heat_of_formation_kjmol": ... }'

4.2 Keshavarz — [Spec-ID: PRED.KESH]

Portée

Corrélations étendues (CHNOFCl, B–N).

Domaines de validité stricts.

AC

 Erreur claire si hors domaine.

 Comparaison sur cartouche test.

Test

pytest -k keshavarz

4.3 ML VoD/Pcj/ΔHdet — [Spec-ID: PRED.ML]

Portée

Descripteurs RDKit + composition + densité + ΔHf.

Modèles RF/GPR/XGBoost ; registry avec hash dataset.

AC

 LOOCV/hold‑out ; métriques R²/MAE/RMSE stockées.

 model_registry.json immuable (hashs).

Test

pytest -k predictors_ml

5) Spectroscopies & analyse (héritage IAM 1.0)
5.1 IR (MVP)

Portée

Backend : POST /api/calc/ir → lance freq (xtb/psi4), parse intensités → CSV/JSON.

Frontend : onglet IR avec graphique interactif (mode barre/ligne).

AC

 CH3NO2 → spectre affiché ; unités (cm⁻¹) ; export CSV.

 Erreurs compréhensibles (pas de fréquences).

Test

curl -X POST /api/calc/ir -d '{ "engine":"psi4", "xyz":"..."}'

5.2 Raman, UV‑Vis, NMR

Portée

Même pipeline/contrat que IR; composant de graphe factorisé.

AC

 4 onglets valides ; axes légendés ; export CSV.

 Données persistées sous Exports/.

Test

pytest -k spectra

6) Visualisation, orbitales & 3D
6.1 Orbitals (cubes)

Portée

Psi4: génération .cube par MO index.

UI : overlay HOMO/LUMO avec 3Dmol.js.

AC

 POST /api/calc/mo_cube → cube_path + viewer URL.

 Lazy‑load + contrôle iso‑value.

Test

pytest -k mo_cubes

6.2 Convert → PDB — [convert/to-pdb]

Portée

RDKit : SMILES/MOL→3D→PDB ; ETKDG + UFF/MMFF.

Retour : {pdb, formula, n_atoms}.

AC

 PDB parseable ; visualisable côté UI.

 Unités/infos moléculaires cohérentes.

Test

curl -X POST /api/convert/to-pdb -d '{"smiles":"CCO"}'

7) Performance CJ (Cantera) & Benchmarks
7.1 Cantera CJ — [predict/cj/v1]

Portée

nox/runners/cantera_cj.py : Pcj, Tcj, composition produits (CHNO + métaux/BN si support).

UI : onglet CJ (tableaux + mini graphiques).

AC

 Endpoint renvoie Pcj, Tcj, fractions produits.

 Export CSV + PNG.

Test

curl -X POST /api/v1/predict/cj -d '{ "formula":"CH3NO2", "rho":..., "Q":... }'

7.2 Benchmarks consolidés

Portée

Outil agrégateur KJ/Emp/ML/CJ avec R², MAE, RMSE ; rapports.

AC

 Rapport .md/.csv généré + graphiques.

 Cartouches tests incluses.

Test

pytest -k benchmark_tool

8) Stabilité & similarité (extensions IAM 1.0)
8.1 Stability Predictor

Portée

SMARTS RDKit : groupes sensibles ; score 0–10 ; conseils.

AC

 Jeu de SMILES de régression → scores cohérents.

 Carte “Stabilité” dans Résumé.

Test

pytest -k stability

8.2 Similarité moléculaire

Portée

Morgan FP + Tanimoto, top‑k voisins depuis IAM_Knowledge.

AC

 Retour voisins + distances ; exemples connus.

Test

pytest -k similarity

9) Ketcher & UI Shell
9.1 Ketcher Integration — [Spec-ID: UI.KETCH] 🔑

Portée

Bundle Ketcher vendored IAM_GUI/static/ketcher/, bridge.js postMessage.

Endpoints utilitaires /ketcher/to-smiles, /ketcher/to-xyz.

AC

 Smoke: /static/ketcher/index.html → 200 ; postMessage OK.

 Sanitation MOL → erreurs 400 cohérentes.

Test

make ketcher.smoke (ou équivalent).

9.2 UI Shell — [Spec-ID: UI.SHELL]

Portée

Tabs: Summary, Output, Orbitals, Spectra, Performance (+ Logs).

Affichage correlation_id dans erreurs ; copy‑to‑clipboard.

AC

 Données peuplées depuis schéma unifié.

 Lazy‑load cubes & spectres.

Test

pytest -k ui_contracts

9.3 Spectra Viewers — [Spec-ID: UI.SPEC]

Portée

IR/Raman/UV‑Vis/NMR -> graphiques + CSV export.

AC

 Axes/Unités corrects ; export CSV.

Test

pytest -k spectra_viewers

10) Agent & automatisation
10.1 IAM_Agent (opt‑in)

Portée

scripts/agent_watch.py : surveille ToAnalyze/ → lance jobs → écrit résultats.

Endpoints agent: /run_shell, /write_file, /log_feedback, /run_backup.

AC

 Boutons GUI déclenchent l’agent ; logs temps réel.

 ZIP de backup listé dans UI.

Test

python scripts/agent_watch.py --dry-run ; pytest -k agent

11) Sécurité, quotas & observabilité

Portée

CORS maîtrisé ; taille max fichiers ; quotas jobs par user.

Compteurs Prometheus (jobs_started/failed, durations).

AC

 429 au‑delà quotas ; 413 si fichier trop gros.

 Métriques exposées (si activé).

Test

pytest -k quotas ; curl /metrics (si présent).

12) Ordre d’implémentation (phases livrables)

M1 – Fondations 🔑
IO.1, IO.2, PIPE.SMILES, XTB (stub + flag), endpoints de conversion, Ketcher bundle (smoke).
Valider: tests de schéma, convert, ketcher smoke.

M2 – XTB réel + OPT→PROP 🔑
XTB réel (flags, timeout), cache OPT, dipôle/gap, artefacts.
Valider: IAM_ENABLE_XTB=1 pytest -k xtb_real.

M3 – Psi4 + Spectres (IR MVP) 🔑
Psi4 SP/OPT/FREQ + IR JSON/CSV ; viewer IR.
Valider: test CH3NO2 IR sur UI + export CSV.

M4 – Prédicteurs
KJ → Keshavarz → ML (registry).
Valider: unités normalisées, métriques ML.

M5 – Orbitals & PDB
MO cubes + convert→PDB.
Valider: 3Dmol overlay HOMO/LUMO + PDB viewer.

M6 – CJ (Cantera) & Benchmarks
Endpoint CJ + rapports consolidés.
Valider: cartouche test Pcj/Tcj.

M7 – Stabilité & Similarité
SMARTS + Morgan/Tanimoto.
Valider: jeux de régression.

M8 – Agent & Historique
Agent, export ZIP, navigateur de résultats.
Valider: bouton GUI → ZIP généré.

13) Tests & CI

Marqueurs pytest: @spec("IO.1"), @spec("ENG.XTB"), … pour lier AC ↔ tests.

CI: ruff + mypy + pytest + “schema conformance”.

External tests: marqués et skippés par défaut (activer via env).
14) Annexes (exemples rapides)

Envelope erreur (standard)
{"code":400,"message":"Invalid input","details":{"field":"smiles"},"correlation_id":"uuid"}
Run XTB (stub/réel selon flag)
curl -s -X POST http://localhost:8000/api/v1/run_xtb \
  -H 'Content-Type: application/json' \
  -d '{"input_type":"smiles","input_data":"CCO","charge":0,"multiplicity":1,"options":{"opt":true}}'
  
IR (Psi4)
IAM_ENABLE_PSI4=1 curl -s -X POST http://localhost:8000/api/calc/ir -d '{"engine":"psi4","xyz":"..."}'
