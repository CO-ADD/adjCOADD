import React, { useState } from 'react'
import { createRoot } from 'react-dom/client'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'
import { Editor } from 'ketcher-react'
import 'ketcher-react/dist/index.css'

// Standalone mode: the indigo struct service is bundled with Ketcher, so no
// remote struct-service URL is required.
const structServiceProvider = new StandaloneStructServiceProvider()

function KetcherEditor({ initialMol }) {
  const [ketcher, setKetcher] = useState(null)

  function handleOnInit(k) {
    setKetcher(k)
    // Expose the instance so page scripts can call getMolfile()/setMolecule().
    window.ketcher = k
    if (initialMol) {
      k.setMolecule(initialMol.replace(/\\n/g, '\n'))
    }
  }

  function updateSmiles() {
    ketcher.getMolfile('v2000').then((molfile) => {
      const out = document.getElementById('id_smol')
      if (out) out.value = molfile
    })
  }

  return (
    <div>
      <div>
        <button id="get_new_mol" type="button" onClick={updateSmiles}>
          Update Smiles
        </button>
      </div>
      <Editor
        staticResourcesUrl=""
        structServiceProvider={structServiceProvider}
        onInit={handleOnInit}
      />
    </div>
  )
}

// Same contract as the old static/js/ketcher bundle:
//   - a hidden <input id="mol"> carries the initial molfile
//   - <div id="root"> is the mount point
//   - the current molfile is written to <input id="id_smol"> on "Update Smiles"
//   - window.ketcher is exposed for external Check/Save handlers
//
// Mounting is explicit (window.initKetcher()) rather than automatic, so the
// bundle can be loaded once per page and re-mounted whenever the modal markup
// is (re)inserted into the DOM.
let currentRoot = null

function initKetcher() {
  const rootEl = document.getElementById('root')
  if (!rootEl) return false

  if (currentRoot) {
    currentRoot.unmount()
    currentRoot = null
  }

  const molInput = document.getElementById('mol')
  const initialMol = molInput ? molInput.value : ''
  currentRoot = createRoot(rootEl)
  currentRoot.render(<KetcherEditor initialMol={initialMol} />)
  return true
}

window.initKetcher = initKetcher
