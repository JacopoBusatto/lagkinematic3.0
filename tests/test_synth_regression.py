from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import yaml
from click.testing import CliRunner
from pandas.testing import assert_frame_equal

from lagkinematic.cli import main as lagk_main


# Percorso del file di riferimento (da tenere nel repo)
REF_CSV = Path(__file__).parent / "data" / "reference_synth_steps.csv"


def _prepare_synth_config(tmp_path: Path) -> Path:
    """
    Copia examples/config_synth.yml in tmp_path e
    sovrascrive solo run.output.dir → tmp_path/OUT_SYNTH_TEST.
    """
    repo_root = Path(__file__).resolve().parents[1]
    src_cfg = repo_root / "examples" / "config_synth.yml"
    if not src_cfg.exists():
        raise FileNotFoundError(f"Config synth non trovato: {src_cfg}")

    with src_cfg.open("r") as f:
        cfg = yaml.safe_load(f)

    outdir = tmp_path / "OUT_SYNTH_TEST"
    outdir.mkdir(parents=True, exist_ok=True)
    cfg["run"]["output"]["dir"] = str(outdir)

    out_cfg = tmp_path / "config_synth_test.yml"
    with out_cfg.open("w") as f:
        yaml.safe_dump(cfg, f)

    return out_cfg


def _load_steps_parquet(outdir: Path) -> pd.DataFrame:
    """
    Legge tutti i Parquet compattati degli steps per tutti i rank
    da outdir/steps/rankXXXX/steps_rankXXXX.parquet
    e li concatena in un singolo DataFrame ordinato per (id, time).
    """
    steps_root = outdir / "steps"
    if not steps_root.exists():
        raise FileNotFoundError(f"Nessuna directory steps trovata in {steps_root}")

    dfs: list[pd.DataFrame] = []

    for rank_dir in sorted(steps_root.glob("rank*")):
        if not rank_dir.is_dir():
            continue

        # Nome generato dal compattatore: steps_rank0000.parquet, ecc.
        parquet_files = sorted(rank_dir.glob("steps_*.parquet"))
        if not parquet_files:
            continue

        table = pq.read_table(parquet_files[0])
        df_rank = table.to_pandas()
        dfs.append(df_rank)

    if not dfs:
        raise RuntimeError(f"Nessun Parquet di steps trovato in {steps_root}")

    df_out = pd.concat(dfs, ignore_index=True)
    df_out = df_out.sort_values(["id", "time"]).reset_index(drop=True)
    return df_out


def _load_reference_df() -> pd.DataFrame:
    """
    Carica il CSV di riferimento per il caso synth,
    ordinato per (id, time).
    """
    if not REF_CSV.exists():
        raise FileNotFoundError(
            f"File di riferimento non trovato: {REF_CSV}\n"
            "Generalo eseguendo lagk3 con config_synth.yml e salvando il Parquet in CSV."
        )

    df_ref = pd.read_csv(REF_CSV)
    df_ref = df_ref.sort_values(["id", "time"]).reset_index(drop=True)
    return df_ref


def test_synth_regression_against_reference(tmp_path):
    """
    Test di regressione per il caso sintetico:

    - copia examples/config_synth.yml in tmp_path e punta l'output a OUT_SYNTH_TEST
    - esegue lagk3
    - legge il Parquet degli step da OUT_SYNTH_TEST/steps
    - confronta il risultato con il CSV di riferimento

    Se in futuro cambiamo integratore/sampling, questo test segnala subito
    eventuali differenze nel comportamento (almeno per il caso sintetico).
    """
    # 1) Prepara config di test
    cfg_path = _prepare_synth_config(tmp_path)

    # 2) Lancia lagk3 tramite CLI runner (senza toccare il filesystem del repo)
    runner = CliRunner()
    result = runner.invoke(lagk_main, ["--config", str(cfg_path)])

    # Se fallisce, stampa l'output CLI nel messaggio di assert
    assert result.exit_code == 0, (
        f"CLI fallita con codice {result.exit_code}:\n{result.output}"
    )

    # 3) Carica output steps
    outdir = tmp_path / "OUT_SYNTH_TEST"
    df_out = _load_steps_parquet(outdir)

    # 4) Carica riferimento
    df_ref = _load_reference_df()

    # 5) Consideriamo solo le colonne comuni, nel caso in futuro aggiungiamo colonne extra
    common_cols = [c for c in df_ref.columns if c in df_out.columns]
    df_out = df_out[common_cols].copy()
    df_ref = df_ref[common_cols].copy()

    # 6) Confronto con tolleranza numerica
    assert_frame_equal(
        df_out,
        df_ref,
        check_dtype=False,
        rtol=1e-7,
        atol=1e-9,
        obj="Synthetic trajectory regression",
    )
