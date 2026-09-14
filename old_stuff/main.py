import argparse

import sys
import os

# functions_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'functions'))
# classes_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'classes'))

# if functions_path not in sys.path:
#     sys.path.append(functions_path)

# if classes_path not in sys.path:
#     sys.path.append(classes_path)

# 3. Import your function directly from the file name
from old_stuff.run_config import RunConfig
from old_stuff.temporal_protocol import TemporalProtocolConfig
from old_stuff.evidence_policy import EvidencePolicy
from old_stuff.selftest import selftest_v24
from old_stuff.run import run_main

def main_v24() -> None:
    """CLI: python depressao_eeg_v24.py --eeg-dir DIR --metadata FILE --output DIR
    --selection-criterion "texto" [--mode research]
    """
    ap = argparse.ArgumentParser(description="MODMA EEG v24: correcoes B1-B11 de auditoria metodologica")
    ap.add_argument("--eeg-dir"); ap.add_argument("--metadata"); ap.add_argument("--output")
    ap.add_argument("--selection-criterion", default="")
    ap.add_argument("--mode", default="integration", choices=["integration", "research"])
    ap.add_argument("--temporal-mode", default="derived", choices=["derived", "fixed"])
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        ok = selftest_v24(True)
        raise SystemExit(0 if all(ok.values()) else 1)
    missing = [k for k in ("eeg_dir", "metadata", "output") if not getattr(args, k)]
    if missing:
        raise SystemExit("Faltam argumentos obrigatorios: %s" % missing)
    cfg = RunConfig(temporal=TemporalProtocolConfig(mode=args.temporal_mode))
    result = run_main(args.eeg_dir, args.metadata, args.output, args.selection_criterion,
                     EvidencePolicy(mode=args.mode), cfg)
    print(json.dumps(result["summary"], indent=2, ensure_ascii=False, default=_json_default)[:4000])


if __name__ == "__main__":
    main_v24()
