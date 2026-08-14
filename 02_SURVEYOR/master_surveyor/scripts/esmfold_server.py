#!/usr/bin/env python
"""ESMFold warm server -- loads the model ONCE and keeps it resident on GPU,
instead of esmfold_predict.py's per-call subprocess (which reloads ~2GB of
weights from disk every time, measured at ~45s of pure overhead per call).
Worth running before any batch that folds more than a handful of sequences.

MUST be invoked with the dedicated `esmfold` conda env's interpreter (same
reason as esmfold_predict.py -- the torch/transformers stack is isolated
there). m2_structures.py's _run_esmfold() calls this over HTTP when
reachable, falling back to the old subprocess-per-call script otherwise, so
forgetting to start this server degrades gracefully rather than failing.

Usage:
    esmfold_server.py [--port 5058]
"""

from __future__ import annotations

import sys
from pathlib import Path

MODEL_NAME = "facebook/esmfold_v1"

_model = None


def _load_model():
    global _model
    if _model is not None:
        return _model
    import torch
    from transformers import EsmForProteinFolding
    print("[esmfold_server] loading model (one-time)...", file=sys.stderr, flush=True)
    _model = EsmForProteinFolding.from_pretrained(MODEL_NAME, low_cpu_mem_usage=True).cuda()
    # Cast just the ESM-2 language-model backbone to half precision (HuggingFace's
    # own documented pattern for facebook/esmfold_v1) -- the backbone alone was
    # eating ~21GB of this box's 24GB in fp32, leaving almost no headroom for
    # activation memory on long sequences (verified: OOM'd on a 2288-residue
    # sequence trying to allocate just 2.5GB more). The folding trunk/structure
    # module stays in fp32 for numerical stability; only the LM half is cast.
    _model.esm = _model.esm.half()
    print("[esmfold_server] model loaded, ready", file=sys.stderr, flush=True)
    return _model


def create_app():
    from flask import Flask, jsonify, request

    app = Flask(__name__)

    @app.route("/health")
    def health():
        return jsonify({"ok": True, "model_loaded": _model is not None})

    @app.route("/fold", methods=["POST"])
    def fold():
        import torch

        body = request.get_json(force=True)
        sequence = body["sequence"]
        out_pdb = Path(body["out_pdb"])
        chunk_size = body.get("chunk_size")

        try:
            model = _load_model()
            model.trunk.set_chunk_size(chunk_size)
            with torch.no_grad():
                pdb_str = model.infer_pdb(sequence)
        except torch.cuda.OutOfMemoryError as exc:
            # Variable-length inputs fragment PyTorch's caching allocator over a
            # long-lived process (verified: total reserved memory grew from ~20GB
            # after one successful fold to ~21.5GB after the next, on a SMALLER
            # sequence) -- empty_cache() releases cached-but-unused blocks back to
            # the driver so the next request isn't starting from a worse position.
            torch.cuda.empty_cache()
            return jsonify({"ok": False, "error": f"CUDA OOM for length {len(sequence)}: {exc}"}), 200
        except Exception as exc:  # noqa: BLE001 -- surfaced to the caller as a normal JSON error, not a 500
            return jsonify({"ok": False, "error": str(exc)}), 200

        torch.cuda.empty_cache()   # same fragmentation mitigation on the success path

        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        out_pdb.write_text(pdb_str)
        return jsonify({"ok": True})

    return app


if __name__ == "__main__":
    port = 5058
    if "--port" in sys.argv:
        port = int(sys.argv[sys.argv.index("--port") + 1])
    _load_model()   # load eagerly at startup, not on first request
    create_app().run(host="127.0.0.1", port=port, threaded=False)   # threaded=False: one GPU, one fold at a time
