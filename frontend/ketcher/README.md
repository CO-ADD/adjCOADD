# Ketcher 3.18.0 build scaffold

Builds a self-contained Ketcher editor bundle into the Django static tree, so the
templates can load it with a plain `<script src>` / `<link>` (no runtime CDN).

## Prerequisites

- Node.js 18+ (the `ketcher-standalone` package declares `node >=24.14.1` for its
  own tooling, but the Vite build itself only needs a recent Node 18/20/22).

## Build

```sh
cd frontend/ketcher
npm install
npm run build
```

Output (written directly into the project static tree):

```
static/js/ketcher3/ketcher.bundle.js    # IIFE bundle (React + ketcher-react + ketcher-standalone)
static/js/ketcher3/styles.css           # extracted ketcher-react styles (Vite lib mode)
static/js/ketcher3/assets/...           # fonts/images/wasm emitted by Ketcher (if any)
```

## Integration contract

The bundle keeps the same contract as the previous `static/js/ketcher/index-bundle.js`:

- `<div id="root">` is the mount point.
- `<input id="mol">` (hidden) carries the initial molfile to load.
- `<input id="id_smol">` (hidden) receives the current molfile when "Update Smiles" is clicked.
- `window.ketcher` is exposed for external Check/Save handlers
  (`window.ketcher.getMolfile('v2000')`, `window.ketcher.setMolecule(...)`).

Mounting is **explicit** — call `window.initKetcher()` after the `#root`/`#mol`
markup is in the DOM (the bundle does not auto-mount). This lets the bundle be
loaded once per page and re-mounted whenever a modal is (re)inserted.

## Templates using it

- `dsample/templates/dsample/coadd/coadd_structure_edit.html`
- `templates/utils/ketcher.html` (used by the Drug structure editor)

## Notes / gotchas

- If Vite emits WASM/worker assets under `static/js/ketcher3/assets/`, and Ketcher
  fails to find its resources, set the `Editor` prop
  `staticResourcesUrl` in `src/main.jsx` to the absolute static URL, e.g.
  `staticResourcesUrl="/static/js/ketcher3/assets/"`.
- A hard browser refresh (Ctrl+F5) is required after swapping the bundle, since the
  old JS/CSS are cached aggressively.
- Do not commit `frontend/ketcher/node_modules/` — add it to `.gitignore` if the
  repo does not already ignore it.
